
/**************************************************START OF FILE*****************************************************/
#ifdef __cplusplus
extern "C"
{
#endif

/*------------------------------------------------------------------------------------------------------------------
Includes
*/
#include "stdint.h"
#include "math.h"//ZykeeAdd
#include "string.h"
#include "stdlib.h"

#include "xALG_FX_Tremolo.h"
#include "xALG_WSOLA_pitch_shift2.h"//ZykeeAdd
#include "fftw3.h"

float parm[10] = {50,50,50,50,50,50,50,50,50,50};


/*------------------------------------------------------------------------------------------------------------------
Code
*/
void Input_BufferChange(float* Input_r , float* Input_l , float* Output , uint32_t Len);
void Output_BufferChange(float* Output_r , float* Output_l , float* Input , uint32_t Len);




/*
*********************************************************************************************************************
@ Brief  : 

@ Param  : NONE

@ Return : NONE

@ Author : YWJ(QQ:872180981)

@ Data   : 2022-07-18 14:42
*********************************************************************************************************************
*/
struct ZykeeFourier ZyF;
void xALG_Basic_Init(void)
{
	xALG_FxTremolo_Init();
    double Pitch = -12;
    ZykeePitchShiftInit(&ZyF,Pitch);
}

/*
*********************************************************************************************************************
@ Brief  : 

@ Param  : NONE

@ Return : NONE

@ Author : YWJ(QQ:872180981)

@ Data   : 2021-10-04 17:59
*********************************************************************************************************************
*/
static float TempBuffer[8 * 1024];
float TransformBuffer[8 * 1024];

void xALG_Basic_Deal(float* InputR , float* InputL , float* Output_R , float* Output_L , uint32_t Len)
{
    float* p = TempBuffer;
    float* q = TransformBuffer;
    Input_BufferChange((float*)InputR, (float*)InputL, p, Len);
    ZykeePitchShift(ZyF.rate, Len, win_len, 4, ZyF.sampleRate,p, q, &ZyF);
	Output_BufferChange(Output_R, Output_L, q, Len);
}


/*
*********************************************************************************************************************
@ Brief  : 左右声道交替

@ Param  : NONE

@ Return : NONE

@ Author : YWJ(QQ:872180981)

@ Data   : 2021-07-09 17:30
*********************************************************************************************************************
*/
void Input_BufferChange(float* Input_r , float* Input_l , float* Output , uint32_t Len)
{
	for (uint32_t i = 0; i < Len; i++) Output[i] = Input_r[i];//ZykeeChanged
}
/*
*********************************************************************************************************************
@ Brief  : 左右声道分开

@ Param  : NONE

@ Return : NONE

@ Author : YWJ(QQ:872180981)

@ Data   : 2021-07-09 17:30
*********************************************************************************************************************
*/
void Output_BufferChange(float* Output_r , float* Output_l , float* Input , uint32_t Len)
{
     for (uint32_t i = 0; i < Len; i++)
    {
         Output_r[i] = Input[i];
         Output_l[i] = Input[i];
    }//ZykeeChanged
}

 
#ifdef __cplusplus
}
#endif

/**************************************************END OF FILE*******************************************************/
