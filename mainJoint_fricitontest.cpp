#include <iostream>
#include <fstream>
#include <string>
#include <vector> 
#include <stdio.h>  
#include <Eigen/Dense>
#include <chrono>
#include <thread>
#include <pthread.h>
#include <sys/time.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>

#ifndef PI
    #define PI 3.141592654
#endif // PI



double frictioncompensation_stribeck(double a, double b, double s,double alpha, double v, double qdot, double ctscale, double gear_ratio);



int main(){
    
    // //
    // TAICHI::twoOrderFilter cartPoleFilter[1];
    // set motor info
    int motorNum = 6;
    Eigen::VectorXd motorIp = Eigen::VectorXd::Zero(motorNum);
    motorIp << 107,25,21,110,109,22;
    // motorIp << 111,106,102,19,103,104;
    // friction test
    // double vmax = 6.0; // rad/s
    Eigen::VectorXd vmax_ = Eigen::VectorXd::Zero(6);
    Eigen::VectorXd vmax = Eigen::VectorXd::Zero(motorNum);
    vmax_ << 5., 5., 5., 5., 5., 5.;
    for(int i = 0;i<motorNum;i++)
    {
        vmax[i] = vmax_[i];
    }    
    double Tduration = 1; // s
    // double vdelt = 0.05; // rad/s
    Eigen::VectorXd vdelt_ = Eigen::VectorXd::Zero(6);
    Eigen::VectorXd vdelt = Eigen::VectorXd::Zero(motorNum);
    vdelt_ << 0.05, 0.05, 0.05, 0.05, 0.05, 0.05;
    for(int i = 0;i<motorNum;i++)
    {
        vdelt[i] = vdelt_[i];
    }   
    // double Ttotal = 0;
    Eigen::VectorXd Ttotal = Eigen::VectorXd::Zero(motorNum);
    Eigen::Matrix<int,Eigen::Dynamic,1> N = Eigen::Matrix<int,Eigen::Dynamic,1>::Zero(motorNum);
    // int N = floor(vmax/vdelt)+1;
    for(int i = 0; i< motorNum;i++)
    {
        N(i) = floor(vmax(i)/vdelt(i))+1;
        Ttotal(i) = Tduration*N(i)*2;
    }
    // Ttotal = Tduration*N*2;
    int simCnt = 0;
    double timeSim = 0.0;
    double timeStep = 0.004;
    int simTotalNum = Ttotal.maxCoeff()/timeStep+10;
    // simTotalNum = 2000;
    double PTPtime = 5;
    // double Ttotal = 60+PTPtime;
    // int simTotalNum = Ttotal/timeStep;
    std::cout<<"simTotalNum: "<<simTotalNum<<std::endl;
    //
    int motionStartCnt = 200;
    int cartMotionCycleCnt = 200;//200--1s
    //
    std::ofstream foutData;
    foutData.open("datacollection.txt",std::ios::out);
    Eigen::VectorXd dataL = Eigen::VectorXd::Zero(3+10*motorNum);

    //----------------------friction compensation--------------------------
    // method 1
    Eigen::MatrixXd para_;
    Eigen::MatrixXd para;
    // ip 111; 106; 102; 19; 103; 104; 107; 25; 21;110; 109; 22
    para_.resize(6,5);
    para_<<  
            1.06064804940437, -5.78371747871328, 11.7627981997819, 111.024140652399, 0.00102757671727644, 
            0.588206167201691, -1.91293190536599, 3.99551564845588, 124.367215752539, -0.00126745046176146,
            0.136497360057077, -1.51004884384212, 2.97111846967942, 35.5769408790899, -0.00111796156728834,
            0.115512410406005, -1.44541794715704, 2.88933111842253, 55.0519847711033, 0.00201362083280513,
            0.216028669928089, -0.906647952260252, 1.89722500778413, 163.227703834085, -0.000262750109186243, 
            0.222786931108455, -0.821761343408838, 1.67176960838594, 87.5358021436918, -5.41660748130355e-05;

            // 0.915262244625368, -5.07648023045021, 10.1374785845289, 137.244164571224, -0.000181470759139340,
            // 0.670457054616744, -1.60338091854170, 3.18061859992825, 71.5111423971947, 0.00156874260830262,
            // 0.132759613590217, -1.53435648835733, 3.10637543035236, 45.3287129894102, 0.00315225307964974,
            // 0.183412133803200, -1.78383773366129, 3.53933989002770, 59.2279900303118, 0.00260733644750306,
            // 0.187301718720714, 0.753668957112634, -1.43505864146487, -119.415622205308, -0.000500017648723363,
            // 0.182227777505046, -0.700032460956267, 1.45590176754402, 102.920723350609, -0.00218173462005282;
    para.resize(motorNum,5);
    for(int i = 0; i < motorNum;i++)
    {
        para.row(i) = para_.row(i);
    }
 
    //----------------------friction compensation--------------------------
    //----------------------identifacation-------------------------------
    Eigen::VectorXd qCmd = Eigen::VectorXd::Zero(motorNum);
    Eigen::VectorXd qDotCmd = Eigen::VectorXd::Zero(motorNum);
    Eigen::VectorXd currentCmd = Eigen::VectorXd::Zero(motorNum);
    Eigen::VectorXd frictioncurrentCmd = Eigen::VectorXd::Zero(motorNum);
    //   
    while (simCnt<simTotalNum)
    {
        
        
        //---------------------high-level control-----------------------
        //
        // friction identification
        int count = floor(timeSim/Tduration);
        for(int i = 0;i <motorNum;i++)
        {
            if(count < N(i))
            {
                qDotCmd[i] = vdelt(i)*count;
                qCmd[i] += qDotCmd[i]*timeStep; 
            }else if(count < 2*N(i))
            {
                qDotCmd[i] = vmax(i) - (count - N(i))*vdelt(i);
                qCmd[i] += qDotCmd[i]*timeStep;
            }else{
                std::cout<<"stop! "<<std::endl;
                qDotCmd[i] = 0.;
                qCmd[i] = qCmd[i];
            }
        }
        
        // --------------------friction compensation-------------------
        ////method 1
        for (int i = 0; i < motorNum; i++){
            // if(abs(qDotEst_kalman[i])>0.1)
            // {
                frictioncurrentCmd[i] = frictioncompensation_stribeck(para(i,0),para(i,1),para(i,2),para(i,3),para(i,4),qDotCmd[i],c_t_scale_[i],gear_[i]);       
            // }else{
            //     frictioncurrentCmd[i] = 0.0;
            // }
            }
        
        for (int i = 0; i < motorNum; i++)
        {
            currentCmd[i] = frictioncurrentCmd[i];
        }
        //---------------------send Cmd--------------------

        
        
        simCnt += 1;
        timeSim =  simCnt*timeStep;
 
        
        
    }
   
    return 0;
}


double frictioncompensation_stribeck(double a, double b, double s,double alpha, double v, double qdot,double ctscale,double gear_ratio)
{
    return (a*qdot + b + s/(1 + exp(-alpha*(qdot + v))))/(ctscale*gear_ratio);
}
