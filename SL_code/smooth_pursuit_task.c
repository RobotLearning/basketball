/*============================================================================
==============================================================================
                      
                              smooth_pursuit_task.c
 
==============================================================================
Remarks:

      follow the first detected blob using neck and eye degrees of freedom

============================================================================*/

// system headers
#include "SL_system_headers.h"

/* SL includes */
#include "SL.h"
#include "SL_user.h"
#include "SL_tasks.h"
#include "SL_task_servo.h"
#include "SL_kinematics.h"
#include "SL_dynamics.h"
#include "SL_collect_data.h"
#include "SL_shared_memory.h"
#include "SL_man.h"
#include "SL_userGraphics.h"

/* defines */

/* local variables */
static int servo_rate = SERVO_BASE_RATE / TASK_SERVO_RATIO;
static int vision_rate, vision_rate_initial = 30;
static char first_vision_time = TRUE;
static double smooth_pursuit_gains_initial[7] = {0.5, 0.1, 1.0, 0.75, 1.0, 0.75, 1.0};
static double blob_change_gains_initial[7] = {0.00025, 0.00015, 0.001, 0.0009, 0.0005, 0.0009, 0.0005};
static double smooth_pursuit_gains[7], blob_change_gains[7];
static char *joint_name[] = {"B_HN", "B_HT", "B_HR", "R_EP", "R_ET", "L_EP", "L_ET"};
static int ball_display = 1;

/* global functions */
void init_smooth_pursuit(void);
int run_smooth_pursuit_task(void);

/* local functions */
static int  init_smooth_pursuit_task(void);
static int  change_smooth_pursuit_task(void);
static double pan_angle(double pixel, double width_pixel);
static double tilt_angle(double pixel, double width_pixel);
static void simulate_ball(int count, float ball[]);

/*****************************************************************************
******************************************************************************
Function Name	: add_smooth_pursuit_task
Date		: Aug. 2010
Remarks:

adds the task to the task menu

******************************************************************************
Paramters:  (i/o = input/output)

none

*****************************************************************************/
void
add_smooth_pursuit_task( void )
{
  int i, j;
  
  addTask("Smooth Pursuit Task", init_smooth_pursuit_task, 
	  run_smooth_pursuit_task, change_smooth_pursuit_task);

}    

void
init_smooth_pursuit( void )
{
int i;
  for (i = 0; i < 7; i++)
  {
    smooth_pursuit_gains[i] = smooth_pursuit_gains_initial[i];
    blob_change_gains[i] = blob_change_gains_initial[i];
  }
  vision_rate = vision_rate_initial;
  first_vision_time = TRUE;
}

/*****************************************************************************
******************************************************************************
  Function Name	: init_smooth_pursuit_task
  Date		: Aug. 2010

  Remarks:

  initialization for task

******************************************************************************
  Paramters:  (i/o = input/output)

       none

 *****************************************************************************/
static int 
init_smooth_pursuit_task(void)
{
int ans;

  init_smooth_pursuit();

  get_int("Display 3-D ball position?", ball_display, &ball_display);

  // ready to go
  ans = 0;
  while (ans != 1)
  {
    if (!get_int("Enter 1 to start or 'q' to abort ...", ans, &ans))
    {
      return FALSE;
    }
  }

  return TRUE;
}

/*****************************************************************************
 ******************************************************************************
 Function Name	: run_smooth_pursuit_task
 Date		: Aug. 2010
 
 Remarks:
 
 task control loop
 
 ******************************************************************************
 Paramters:  (i/o = input/output)
 
 none
 
 *****************************************************************************/
int
run_smooth_pursuit_task()
{
// Joint limits
const double min_max_angles[7][2] = {{-0.4f, 0.44f}, {-0.22f, 0.22f},
                                     {-0.79f, 0.79f}, {-0.65f, 0.65f},
                                     {-0.65f, 0.65f}, {-0.65f, 0.65f},
                                     {-0.65f, 0.65f}};
// Velocity limits
const double max_step_eyes = 0.005f, max_step_head = 0.0025f;

int normal_time_steps = servo_rate / vision_rate;
double factor = ((float) servo_rate) / ((float) vision_rate);
double width = 640, height = 480;
int last_frame_counter = -1;

double update[7][2], tmp_update[7];

static int time_steps = 0, extra_time_steps = 0;
static float extra_update[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
             time_factor = 1.0f;
static float previous_blob_left[2], previous_blob_right[2],
             blob_change_left[2], blob_change_right[2];
static float new_state[7], previous_state[7], current_update[7];

int i, j;

float ball[N_CART];

  if (ball_display)
  {
    if (real_robot_flag)
    {
      if (blobs[1].status == 1)
      // Transform result into world coordinates
        for (i = 1; i <= 3; i++)
        {
          ball[i-1] = Alink[BASE][i][4];
          for (j = 1; j <= 3; j++)
            ball[i-1] += Alink[BASE][i][j] * blobs[1].blob.x[j];
        }
    }
    else
    {
    // Simulate ball motion
    static unsigned int count = 0;
      simulate_ball(count, ball);
      count++;
      normal_time_steps = 1; // factor = 8;
    }
    if (blobs[1].status == 1)
      sendUserGraphics((char *) "ball", ball, sizeof(ball));
  }

  if (((raw_blobs2D[1][1].status == 1 || raw_blobs2D[1][2].status == 1) &&
      !first_vision_time) ||
      (raw_blobs2D[1][1].status == 1 && raw_blobs2D[1][2].status == 1))
  {
  // New vision blob has arrived.
    if (first_vision_time)
    {
      first_vision_time = FALSE;
      for (i = 0; i < 2; i++)
      {
        blob_change_left[i] = blob_change_right[i] = 0;
        previous_blob_left[i] = raw_blobs2D[1][1].x[_X_+i];
        previous_blob_right[i] = raw_blobs2D[1][2].x[_X_+i];
      }
    }
    else
    {
      for (i = 0; i < 2; i++)
      {
        blob_change_left[i] = raw_blobs2D[1][1].x[_X_+i] - previous_blob_left[i];
        blob_change_right[i] = raw_blobs2D[1][2].x[_X_+i] - previous_blob_right[i];
        previous_blob_left[i] = raw_blobs2D[1][1].x[_X_+i];
        previous_blob_right[i] = raw_blobs2D[1][2].x[_X_+i];
      }
    }

    time_steps = extra_time_steps = 0;
    time_factor = 1.0f;
    for (i = 0; i < 7; i++)
      extra_update[i] = 0.0f;

    for (i = B_HN; i <= B_HR; i++)
      previous_state[i-B_HN] = joint_des_state[i].th;
      // previous_state[i-B_HN] = joint_state[i].th;
    for (i = R_EP; i <= L_ET; i++)
      previous_state[i-B_HN] = joint_des_state[i].th;
      // previous_state[i-B_HN] = joint_state[i].th;

    /*****************************************************************************

     Calculate the updates for head and eye joints

    *****************************************************************************/
    if (raw_blobs2D[1][1].status == 1)
    {
      update[0][0] = smooth_pursuit_gains[0] *
                     tilt_angle(raw_blobs2D[1][1].x[_Y_], height) +
                     blob_change_gains[0] * blob_change_left[1];
      update[6][0] = smooth_pursuit_gains[6] *
                     tilt_angle(raw_blobs2D[1][1].x[_Y_], height) +
                     blob_change_gains[6] * blob_change_left[1];

      update[5][0] = (smooth_pursuit_gains[5] *
                      pan_angle(raw_blobs2D[1][1].x[_X_], width)) -
                     blob_change_gains[5] * blob_change_left[0];
      // Switched 2 and 1 -- seems different in SL
      update[1][0] = (smooth_pursuit_gains[1] *
                      pan_angle(raw_blobs2D[1][1].x[_X_], width)) -
                     blob_change_gains[1] * blob_change_left[0];
      update[2][0] = smooth_pursuit_gains[2] *
                     pan_angle(raw_blobs2D[1][1].x[_X_], width) +
                     blob_change_gains[2] * blob_change_left[0];
    }

    if (raw_blobs2D[1][2].status == 1)
    {
      update[0][1] = smooth_pursuit_gains[0] *
                     tilt_angle(raw_blobs2D[1][2].x[_Y_], height) +
                     blob_change_gains[0] * blob_change_right[1];
      update[4][1] = smooth_pursuit_gains[4] *
                     tilt_angle(raw_blobs2D[1][2].x[_Y_], height) +
                     blob_change_gains[4] * blob_change_right[1];

      update[3][1] = (smooth_pursuit_gains[3] *
                      pan_angle(raw_blobs2D[1][2].x[_X_], width)) -
                      blob_change_gains[3] * blob_change_right[0];
      // Switched 2 and 1 -- seems different in SL
      update[1][1] = (smooth_pursuit_gains[1] *
                      pan_angle(raw_blobs2D[1][2].x[_X_], width)) -
                     blob_change_gains[1] * blob_change_right[0];
      update[2][1] = smooth_pursuit_gains[2] *
                     pan_angle(raw_blobs2D[1][2].x[_X_], width) +
                     blob_change_gains[2] * blob_change_right[0];
    }

    /*****************************************************************************
       
     Combine the updates for left and right eye
      
    *****************************************************************************/
    if (raw_blobs2D[1][1].status == 1)
    {
      if (raw_blobs2D[1][2].status != 1)
      {
        for (i = 0; i < 7; i++)
          current_update[i] = update[i][0];
        current_update[3] = current_update[5] - previous_state[3] + previous_state[5];
        current_update[4] = current_update[6] - previous_state[4] + previous_state[6];
      }
      else
      {
        for (i = 0; i < 3; i++)
          current_update[i] = (update[i][0] + update[i][1]) / 2.0f;
        current_update[3] = update[3][1];
        current_update[4] = update[4][1];
        current_update[5] = update[5][0];
        current_update[6] = update[6][0];
      }
    }
    else if (raw_blobs2D[1][2].status == 1)
    {
      for (i = 0; i < 7; i++)
        current_update[i] = update[i][1];
      current_update[5] = current_update[3] - previous_state[5] + previous_state[3];
      current_update[6] = current_update[4] - previous_state[6] + previous_state[4];
    }

    current_update[5] = -current_update[5];

    // Take into account that servo rate is higher than vision frame rate
    for (i = 0; i < 7; i++)
      current_update[i] /= factor;

    // Check the velocity limits
    for (i = 0; i < 3; i++)
    {
      if (current_update[i] > max_step_head)
        current_update[i] = max_step_head;
      else if (current_update[i] < -max_step_head)
        current_update[i] = -max_step_head;
    }
    for (i = 3; i < 7; i++)
    {
      if (current_update[i] > max_step_eyes)
        current_update[i] = max_step_eyes;
      else if (current_update[i] < -max_step_eyes)
        current_update[i] = -max_step_eyes;
    }
  }

  // Start stopping the movement if vision frames are lost
  if (time_steps == normal_time_steps)
  {
    extra_time_steps++;
    time_factor /= 2;
    for (i = 0; i < 7; i++)
    {
      extra_update[i] = current_update[i] * time_factor;
      tmp_update[i] = normal_time_steps * current_update[i] +
                      extra_update[i];
    }
  }
  else
  {
    time_steps++;
    for (i = 0; i < 7; i++)
      tmp_update[i] = time_steps * current_update[i];
      // tmp_update[i] = normal_time_steps * current_update[i];
  }

  // Caluclate joint desired state and check the joint limits
  for (i = 0; i < 7; i++)
  {
    new_state[i] = previous_state[i] + tmp_update[i];

    if (new_state[i] < min_max_angles[i][0])
      new_state[i] = min_max_angles[i][0];
    else if (new_state[i] > min_max_angles[i][1])
      new_state[i] = min_max_angles[i][1];
  }

  // Apply the new desired joint configuration
  for (i = R_EP; i <= L_ET; i++)
    joint_des_state[i].th = new_state[i-B_HN];
  for (i = B_HN; i <= B_HR; i++)
    joint_des_state[i].th = new_state[i-B_HN];

  /* Mark that the blob information has been used. Don't use it before the next
     blob arrives. */
  raw_blobs2D[1][1].status = raw_blobs2D[1][2].status = 0;

  if (frame_counter != last_frame_counter)
  {
    last_frame_counter = frame_counter;
    // printf("Status: %d %d\n", count, frame_counter);
  }
    
  return TRUE;
}

/*****************************************************************************
******************************************************************************
  Function Name	: change_sample_task
  Date		: Dec. 1997

  Remarks:

  changes the task parameters

******************************************************************************
  Paramters:  (i/o = input/output)

  none

 *****************************************************************************/
static int 
change_smooth_pursuit_task(void)
{
int i, j1 = 0, ivar;
double dvar;
unsigned char cont = TRUE;
char text[256];

  while (cont)
  {
    j1 = 0;
    printf("\n");
    for (i = B_HN; i <= L_ET; i++)
      printf("%3d %5s gain           %f\n", i - B_HN + 1, joint_name[i-B_HN],
              smooth_pursuit_gains[i-B_HN]);
    printf("\n");
    for (i = B_HN; i <= L_ET; i++)
      printf ("%3d %5s velocity gain  %f\n", i - B_HN + 1 + (L_ET-B_HN+1),
              joint_name[i-B_HN], blob_change_gains[i-B_HN]);
    printf("\n%3d  vision frame rate   %d\n", 2*(L_ET-B_HN+1) + 1, vision_rate);
    printf("%3d  display ball (1/0)  %d\n", 2*(L_ET-B_HN+1) + 2, ball_display);

    get_int("\nEnter parameter index, or <return> to exit", j1, &j1);

    if (j1 >= 1 && j1 <= L_ET-B_HN+1)
    {
      sprintf(text, "Enter new gain for %s", joint_name[j1-1]);
      get_double(text, smooth_pursuit_gains[j1-1], &(smooth_pursuit_gains[j1-1]));
    }
    else if (j1 >= L_ET-B_HN+2 && j1 <= 2*(L_ET-B_HN+1))
    {
      sprintf(text, "Enter new blob chain gain for %s", joint_name[j1-1-(L_ET-B_HN+1)]);
      get_double(text, blob_change_gains[j1-1-(L_ET-B_HN+1)],
                 &(blob_change_gains[j1-1-(L_ET-B_HN+1)]));
    }
    else if (j1 == 2*(L_ET-B_HN+1)+1)
    {
      sprintf(text, "Enter new vision frame rate");
      get_int(text, vision_rate, &(vision_rate));
    }
    else if (j1 == 2*(L_ET-B_HN+1)+2)
    {
      sprintf(text, "Do you want to display ball? (1/0)");
      get_int(text, ball_display, &(ball_display));
      if (!ball_display)
        sendUserGraphics((char *) "clear", NULL, 0);
    }
    else
      cont = FALSE;
  }

  return TRUE;
}

static double pan_angle(double pixel, double width_pixel)
{
  double pixel_offset_horiz = pixel - (width_pixel - 1.0) / 2.0;
  return (0.0005 * pixel_offset_horiz);
}

static double tilt_angle(double pixel, double height_pixel)
{
  double foveation_offset = 0; // 30;
  double pixel_offset_vert = pixel - (height_pixel - 1.0) / 2.0 + foveation_offset;
  return (0.0005 * pixel_offset_vert);
}

/*****************************************************************************
 ******************************************************************************
 Function Name	: init_smooth_pursuit_task
 Date		: Aug. 2010
 
 Remarks:
 
 Simulate 3-D ball motion and project it onto the stereo camera image using
 similar camera parameters as cb-i
 
 ******************************************************************************
 Paramters:  (i/o = input/output)
 
 none
 
 *****************************************************************************/
/*  */ 
static void simulate_ball(int count, float ball[])
{
double ball_eye[N_CART+1], tmp;
int i, j;

  // Ball motion in world coordinates.
  ball[0] = (sin(count/1000.0*3.141592)+1)/2*0.436;
  ball[1] = 1.5; // (sin(count/500.0*3.141592)+1)/2*0.484;
  ball[2] = (sin(count/500.0*3.141592)+1)/2*0.529;
  for (i = 1; i <= 3; i++)
  {
    blobs[1].blob.x[i] = ball[i-1];
  } 
  blobs[1].status = 1;

  // Transform world coordinates to left eye coordinates
  for (i = 1; i <= 3; i++)
  {
    ball_eye[i] = 0;
    for (j = 1; j <= 3; j++)
      ball_eye[i] += Alink[L_EYE_AXIS][j][i] * (blobs[1].blob.x[j] - Alink[L_EYE_AXIS][j][4]);
    ball[i-1] = blobs[1].blob.x[i];
  }

  // Project the ball onto the left eye camera image.
  tmp = ball_eye[1];
  ball_eye[1] = -ball_eye[3];
  ball_eye[3] = -ball_eye[2];
  ball_eye[2] = tmp;
  raw_blobs2D[1][1].status = 1;
  raw_blobs2D[1][1].x[_X_] = 520 * ball_eye[1] / ball_eye[3] + 320;
  raw_blobs2D[1][1].x[_Y_] = 520 * ball_eye[2] / ball_eye[3] + 240;
    
  // Transform world coordinates to right eye coordinates
  for (i = 1; i <= 3; i++)
  {
    ball_eye[i] = 0;
    for (j = 1; j <= 3; j++)
      ball_eye[i] += Alink[R_EYE_AXIS][j][i] * (blobs[1].blob.x[j] - Alink[R_EYE_AXIS][j][4]);
    ball[i-1] = blobs[1].blob.x[i];
  }
    
  // Project the ball onto the right eye camera image.
  tmp = ball_eye[1];
  ball_eye[1] = -ball_eye[3];
  ball_eye[3] = -ball_eye[2];
  ball_eye[2] = tmp;
  raw_blobs2D[1][2].status = 1;
  raw_blobs2D[1][2].x[_X_] = 520 * ball_eye[1] / ball_eye[3] + 320;
  raw_blobs2D[1][2].x[_Y_] = 520 * ball_eye[2] / ball_eye[3] + 240;
}


