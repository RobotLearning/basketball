/*============================================================================
==============================================================================
                      
                              initUserGraphics.c
 
==============================================================================
Remarks:

         Functions needed for user graphics
         simulation

============================================================================*/

#include "SL.h"
#include "SL_user.h"
#include "SL_man.h"

// openGL includes
#include "GL/glut.h"
#include "SL_openGL.h"
#include "SL_userGraphics.h"
#include "math.h"

/**
 * Adds stripes to the raw basketball sphere
 */
void draw_stripes(double b[N_CART+1]) {

	double basketball_radius = b[0];
	int lineAmount = 100; //# of triangles used to draw circle
	GLfloat thickness = 2.0;
	int i;
	double x,y;
	glLineWidth(thickness);
	// draw horizontal stripe
	glBegin(GL_LINE_LOOP);
		for(i = 0; i <= lineAmount;i++) {
			glVertex3f( b[1] + (basketball_radius * cos(i*2*PI/lineAmount)),
			            b[2] + (basketball_radius * sin(i*2*PI/lineAmount)),
				        b[3]);
		}
	glEnd();

	// draw vertical stripe
	glBegin(GL_LINE_LOOP);
		for(i = 0; i <= lineAmount;i++) {
			glVertex3f( b[1],
			            b[2] + (basketball_radius * cos(i*2*PI/lineAmount)),
				        b[3] + (basketball_radius * sin(i*2*PI/lineAmount)));
		}
	glEnd();

	double angle = PI * 30.0/180;
	// draw two stripes below and above horizontal
	double r = basketball_radius * cos(angle);
	glBegin(GL_LINE_LOOP);
		for(i = 0; i <= lineAmount;i++) {
			glVertex3f(b[1] + (r*cos(i*2*PI/lineAmount)),
					b[2] + (r*sin(i*2*PI/lineAmount)), b[3] + basketball_radius*sin(angle));
		}
	glEnd();
	glBegin(GL_LINE_LOOP);
		for(i = 0; i <= lineAmount;i++) {
			glVertex3f(b[1] + (r*cos(i*2*PI/lineAmount)),
					b[2] + (r*sin(i*2*PI/lineAmount)), b[3] - basketball_radius*sin(angle));
		}
	glEnd();
	glLineWidth(1.f);//reset
}

/**
 * Draws the basketball without stripes
 */
void draw_ball(double b[N_CART+1]) {

	double basketball_radius = b[0];
	double basketball_colors[3] = {207.0/256, 83.0/256, 0.0/256};
	GLfloat col[N_CART+1] = {(float)basketball_colors[0],
			                 (float)basketball_colors[1],
							 (float)basketball_colors[2],
							 (float)1.0};

	glPushMatrix();
	glTranslated((GLdouble)b[1],(GLdouble)b[2],(GLdouble)b[3]);
	glColor4fv(col);
	glMaterialfv(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,col);
	glutSolidSphere(basketball_radius,8,8);
	glPopMatrix();
}

/**
 * Draws the ball and puts 1 vertical, 3 horizontal stripes on it
 */
void display_basketball(void *b) {

	double vars[N_CART+1];
	memcpy(&(vars[0]),b,(N_CART+1)*sizeof(double));
	draw_ball(vars);
	draw_stripes(vars);

}

/*
 * Draws the basketball and attaches it to a string
 */
void display_basketball_pendulum(void *data) {

	int i;
	double string_len;
	double angle;
	double ball_radius;
	double vars[6];
	double base_pendulum[N_CART];

	// extract the data
	memcpy(&vars[0],data,6*sizeof(double));
	string_len = vars[0];
	ball_radius = vars[1];
	base_pendulum[0] = vars[2];
	base_pendulum[1] = vars[3];
	base_pendulum[2] = vars[4];
	angle = vars[5];

	double ball_info[N_CART+1];
	ball_info[0] = ball_radius;
	ball_info[1] = base_pendulum[0]; //x is fixed
	ball_info[2] = base_pendulum[1] - (string_len + ball_radius) * sin(angle);
	ball_info[3] = base_pendulum[2] - (string_len + ball_radius) * cos(angle);

	display_basketball(ball_info); // send radius and centre info
	// draw the pendulum string
	glBegin(GL_LINE_LOOP);
	glVertex3f(base_pendulum[0],base_pendulum[1],base_pendulum[2]);
	glVertex3f(base_pendulum[0],base_pendulum[1] - string_len*sin(angle), base_pendulum[2] - string_len*cos(angle));
	glEnd();
}

// local variables

/*****************************************************************************
******************************************************************************
Function Name	: initUserGraphics
Date		: June 1999
   
Remarks:

      allows adding new graphics functions to openGL interface

******************************************************************************
Paramters:  (i/o = input/output)

  none   

*****************************************************************************/
int
initUserGraphics(void)

{
	addToUserGraphics("basketball","Display basketball", &(display_basketball), 4*sizeof(double));
	addToUserGraphics("basketball_pendulum","Basketball attached to string", &(display_basketball_pendulum), 6*sizeof(double));
	return TRUE;

}

