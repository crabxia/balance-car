/***********************************************
公司：轮趣科技（东莞）有限公司
品牌：WHEELTEC
官网：wheeltec.net
淘宝店铺：shop114407458.taobao.com 
速卖通: https://minibalance.aliexpress.com/store/4455017
版本：5.7
修改时间：2021-04-29

 
Brand: WHEELTEC
Website: wheeltec.net
Taobao shop: shop114407458.taobao.com 
Aliexpress: https://minibalance.aliexpress.com/store/4455017
Version:5.7
Update：2021-04-29

All rights reserved
***********************************************/
#ifndef __CONTROL_H
#define __CONTROL_H
#include "sys.h"

extern int Sensor_Left,Sensor_Middle,Sensor_Right,Sensor;
#define PI 3.14159265							//PI圆周率
#define Control_Frequency  200.0	//编码器读取频率
#define Diameter_67  67.0 				//轮子直径67mm 
#define EncoderMultiples   4.0 		//编码器倍频数
#define Encoder_precision  13.0 	//编码器精度 13线
#define Reduction_Ratio  30.0			//减速比30
#define Perimeter  210.4867 			//周长，单位mm

//小车各模式定义
#define Normal_Mode							0
#define Ultrasonic_Avoid_Mode               1
#define Ultrasonic_Follow_Mode              2
#define Lidar_Avoid_Mode					3
#define Lidar_Follow_Mode					4
#define Lidar_Straight_Mode                 5
#define ELE_Line_Patrol_Mode				7
#define CCD_Line_Patrol_Mode				6
#define Move_Specific_Distance_Mode			8 //移动特定距离模式

//避障模式的参数
#define  avoid_Distance 350//避障距离350mm
#define avoid_Angle1 50 //避障的角度，在310~360、0~50°的范围
#define avoid_Angle2 310
#define avoid_speed 30    //避障速度
#define turn_speed 1000;    //避障转向速度

//雷达走直线的参数
#define Initial_speed 30//小车的初始速度大概为200mm每秒
#define Limit_time 600   //限制时间，5ms中断*数值=时间 ，这里就是3s
#define refer_angle1  71 //参照物的角度1
#define refer_angle2  74 //参照物的角度2

//雷达跟随参数
#define Follow_distance 1500  //雷达跟随模式最远距离
#define tracking_speed 40      //给小车一个大概300mm/s的速度


#define Detect_distance 700//检测距离为700mm
#define Barrier_Detected						1
#define No_Barrier								0

#define DIFFERENCE 100

int EXTI15_10_IRQHandler(void);
int Balance(float angle,float gyro);
int Velocity(int encoder_left,int encoder_right);
int Turn(float gyro);
void Set_Pwm(int motor_left,int motor_right);
void Key(void);
void Limit_Pwm(void);
int PWM_Limit(int IN,int max,int min);
u8 Turn_Off(float angle, int voltage);
void Get_Angle(u8 way);
int myabs(int a);
int Pick_Up(float Acceleration,float Angle,int encoder_left,int encoder_right);
int Put_Down(float Angle,int encoder_left,int encoder_right);
void Get_Velocity_Form_Encoder(int encoder_left,int encoder_right);
u8 Choose(void);
void Lidar_Avoid(void);
void Lidar_Straight(void);
void Lidar_Follow(void);
void  Find_CCD_Zhongzhi(void);
u8 Detect_Barrier(void);
void CCD_Mode(void);
int CCD_turn(u8 CCD,float gyro);//转向控制
int ELE_turn(float gyro);//转向控制
void ELE_Mode(void);
void Select_Zhongzhi(void);

int Position_FeedbackControl(float TargetPosition, float CurrentPosition);
int PWM_Restrict(int PWM_P, float TargetVelocity);
int Velocity_FeedbackControl(int PWM_V, int Encoder);
int Line(int CurrentPosition_Left, int CurrentPosition_Right );

#endif
