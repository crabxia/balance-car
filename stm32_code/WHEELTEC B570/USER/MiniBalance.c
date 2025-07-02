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
#include "stm32f10x.h"
#include "sys.h"
#include "control.h"
u8 Way_Angle = 2;																								 // 获取角度的算法，1：四元数  2：卡尔曼  3：互补滤波
u8 Flag_front, Flag_back, Flag_Left, Flag_Right, Flag_velocity = 2;												 // 蓝牙遥控相关的变量
u8 Flag_Stop = 1, Flag_Show = 0;																				 // 电机停止标志位和显示标志位  默认停止 显示打开
int Motor_Left, Motor_Right;																					 // 电机PWM变量 应是Motor的 向Moto致敬
int Temperature;																								 // 温度变量
int Voltage, Middle_angle;																				 // 电池电压采样相关的变量、机械中值
u8 Mode = 0;																									 // 模式选择，默认是普通的控制模式
u8 CCD_Zhongzhi, CCD_Yuzhi;																						 // CCD中值和阈值
u8 LD_Successful_Receive_flag;																					 // 雷达成功接收数据标志位
u16 ADV[128] = {0};																								 // 存放CCD的数据的数组
u16 determine;																									 // 雷达跟随模式的一个标志位
u8 Lidar_Detect = Lidar_Detect_ON;																				 // 电磁巡线模式雷达检测障碍物，默认开启
float  Gyro_Balance, Gyro_Turn;																	 // 平衡倾角 平衡陀螺仪 转向陀螺仪
u32 Distance;																									 // 超声波测距
u8 delay_50, delay_flag, PID_Send;																				 // 延时和调参相关变量
u8 Flag_follow = 0, Flag_avoid = 0;																				 // 超声波跟随、超声波壁障标志位
float Acceleration_Z;																							 // Z轴加速度计
float RC_Velocity, RC_Turn_Velocity, Move_X, Move_Z, PS2_ON_Flag;												 // 遥控控制的速度
float  Velocity_Kp = 16000, Velocity_Ki = 80; // PID参数（放大100倍）
u8 one_lap_data_success_flag = 0;																				 // 雷达数据完成一圈的接收标志位
int lap_count = 0;																								 // 当前雷达这一圈数据有多少个点
int PointDataProcess_count = 0, test_once_flag = 0, Dividing_point = 0;											 // 雷达接收数据点的计算值、接收到一圈数据最后一帧数据的标志位、需要切割数据的数据数



//当前的平衡速度，当前的平衡角角加速度
float Angle_Balance,Gyro_Balance;
int Encoder_Left, Encoder_Right; 
int CurrentPosition_Left=0, CurrentPosition_Right=0,CurrentPosition=0,CurrentPosition_Last=0; // 位置
float TargetPosition=0,TargetVelocity; //实现先快后慢
float Velcity_Kp=80,  Velcity_Ki=0.1,  Velcity_Kd=0; //相关速度PID参数
float Position_Kp=1, Position_Ki=0.01, Position_Kd=0.3; //相关位置PID参数
float Balance_Kp = 325, Balance_Kd = 1.3;
float Turn_Kp = 42, Turn_Kd = 0.6;
float Line_Kp, Line_Kd;

float Acceleration_X;

// float Velcity_Kp=0.4,  Velcity_Ki=0.05,  Velcity_Kd=0; //相关速度PID参数
// float Position_Kp=0.3, Position_Ki=0, Position_Kd= 0; //相关位置PID参数
// float Balance_Kp = 150, Balance_Kd = 0.8;
// float Turn_Kp = 2.5, Turn_Kd = 0;

// float Velcity_Kp=0,  Velcity_Ki=0,  Velcity_Kd=0; //相关速度PID参数
// float Position_Kp=0, Position_Ki=0, Position_Kd=0; //相关位置PID参数
// float Balance_Kp = 25500, Balance_Kd = 135;
// float Turn_Kp = 0, Turn_Kd = 0;
int main(void)
{
	MY_NVIC_PriorityGroupConfig(2); // 设置中断分组
	delay_init();					// 延时函数初始化
	JTAG_Set(JTAG_SWD_DISABLE);		// 关闭JTAG接口
	JTAG_Set(SWD_ENABLE);			// 打开SWD接口 可以利用主板的SWD接口调试
	LED_Init();						// 初始化与 LED 连接的硬件接口
	KEY_Init();						// 按键初始化
	MiniBalance_PWM_Init(7199, 0);	// 初始化PWM 10KHZ与电机硬件接口，用于驱动电机
	uart_init(115200);				// 串口1初始化
	uart3_init(9600);				// 串口3初始化，用于蓝牙模块
	Encoder_Init_TIM3();			// 编码器接口
	Encoder_Init_TIM4();			// 初始化编码器4
	Adc_Init();						// adc初始化
	IIC_Init();						// IIC初始化
	OLED_Init();					// OLED初始化
	//    PS2_Init();                   //ps2手柄初始化,PS2和CCD巡线、ELE巡线不能同时使用
	MPU6050_initialize(); // MPU6050初始化
	DMP_Init();			  // 初始化DMP
	// while(Choose()) { }
	// if(Mode==Ultrasonic_Avoid_Mode||Mode==Ultrasonic_Follow_Mode)
	// 	TIM2_Cap_Init(0XFFFF,72-1); //超声波和雷达、CCD、ELE巡线不能同时使用，使用CCD、ELE功能时，需要拆下超声波模块
	// if(Mode==Lidar_Avoid_Mode||Mode==Lidar_Follow_Mode||Mode==Lidar_Straight_Mode)
	// 	LIDAR_USART_Init();
	// if(Mode==CCD_Line_Patrol_Mode)
	// 	ccd_Init();
	// if(Mode==ELE_Line_Patrol_Mode)
	// 	ele_Init();
	Mode = 8;
	


	MiniBalance_EXTI_Init(); // MPU6050 5ms定时中断初始化，节省定时器资源，减少cpu负担
	while (1)
	{
		if (Flag_Show == 0) // 使用MiniBalance APP和OLED显示屏
		{
			PS2_Read();
			APP_Show();	 // 发送数据给APP
			oled_show(); // 显示屏打开
		}
		else // 使用MiniBalance上位机 上位机使用的时候需要严格的时序，故此时关闭app监控部分和OLED显示屏
		{
			DataScope(); // 开启MiniBalance上位机
		}
		delay_flag = 1;
		delay_50 = 0;
		while (delay_flag)
			; // 示波器需要50ms	高精度延时，delay函数不满足要求，故使用MPU6050中断提供50ms延时
	}
}
