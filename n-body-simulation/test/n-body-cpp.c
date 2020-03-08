/*
 * bulirsch_stoer.cpp
 *
 * Copyright 2011-2013 Mario Mulansky
 * Copyright 2011-2012 Karsten Ahnert
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or
 * copy at http://www.boost.org/LICENSE_1_0.txt)
 */


 /* In order to accomplish this project I took a benchmark program from the above listed boost library.
 This code contain subroutines that were not written by me, however I do understand what they are used
 for and how to use them. The things that I added to this routine are the functions that contain and
obtain the mass for different planets as well as the initial conditions and all of the routines to compute
the derivatives. In order to perform error checking I wrote my own routine that checks the conservation
energy instead of using and built in error monitoring routines.


 RC
 */

#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <omp.h>

#include <boost/array.hpp>
#include <boost/ref.hpp>

#include <boost/numeric/odeint/config.hpp>

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer_dense_out.hpp>

const double newtonsg = 6.647e-11;

using namespace std;
using namespace boost::numeric::odeint;

const int number_bodies=7;  //number of bodies in the system


typedef boost::array< double , (number_bodies*6) > state_type;  //the following are arrays containing
typedef boost::array< double , number_bodies > state_type_1;    //all the different data tupes
typedef boost::array< double , number_bodies > state_type_2;
typedef boost::array< double , number_bodies > state_type_3;
typedef boost::array< double , number_bodies > state_type_4;

/*
 * x' = ( - x*sin t  + 2 tan x ) y
 * with x( pi/6 ) = 2/sqrt(3) the analytic solution is 1/cos t
 */
    //double bot=0;

 double get_denom(const state_type &x , int one, int two )  //routine that computes distance berween bodies
 {
        int jb=6*two;
        int ib=6*one;
        double tmp_denom=pow(x[jb]-x[ib],2)+pow(x[jb+2]-x[ib+2],2)+pow(x[jb+4]-x[ib+4],2);
        return tmp_denom;

 }

 double mass(int one)   //function that provides the mass
 {
        double mass[number_bodies]={0.961*1.989e30,0.006*1.898e27,0.009*1.898e27,0.023*1.898e27,0.025*1.898e27,0.006*1.898e27,0.079*1.898e27};
        return mass[one];
 }

void rhs( const state_type &x , state_type &dxdt , const double t ) //function that solves the derivative
{   //I unfortunately had to do these 1 by 1 as the routine I had to do this process with loops refused to work.

    double denom1=pow(sqrt(get_denom(x,0,1)),3); // sun and earth
    double denom2=pow(sqrt(get_denom(x,0,2)),3); // sun and mars
    double denom3=pow(sqrt(get_denom(x,1,2)),3); // earth and mars
    double denom4=pow(sqrt(get_denom(x,0,3)),3); // sun and jupiter
    double denom5=pow(sqrt(get_denom(x,1,3)),3); // earth and jupiter
    double denom6=pow(sqrt(get_denom(x,2,3)),3); // mars and jupiter
    double denom7=pow(sqrt(get_denom(x,0,4)),3); //sun and saturn
    double denom8=pow(sqrt(get_denom(x,1,4)),3); //earth and saturn
    double denom9=pow(sqrt(get_denom(x,2,4)),3); //mars and saturn
    double denom10=pow(sqrt(get_denom(x,3,4)),3); //jupiter and saturn
    double denom11=pow(sqrt(get_denom(x,0,5)),3); //sun and uranus
    double denom12=pow(sqrt(get_denom(x,1,5)),3); //earth and uranus
    double denom13=pow(sqrt(get_denom(x,2,5)),3); //mars and uranus
    double denom14=pow(sqrt(get_denom(x,3,5)),3); //jupiter and uranus
    double denom15=pow(sqrt(get_denom(x,4,5)),3); //saturn and uranus
    double denom16=pow(sqrt(get_denom(x,0,6)),3); //sun and neptune
    double denom17=pow(sqrt(get_denom(x,1,6)),3); //earth and neptune
    double denom18=pow(sqrt(get_denom(x,2,6)),3); //mars and neptune
    double denom19=pow(sqrt(get_denom(x,3,6)),3); //jupiter and neptune
    double denom20=pow(sqrt(get_denom(x,4,6)),3); //saturn and neptune
    double denom21=pow(sqrt(get_denom(x,5,6)),3);  //uranus and neptune

    double mass1=mass(0);
    double mass2=mass(1);
    double mass3=mass(2);
    double mass4=mass(3);
    double mass5=mass(4);
    double mass6=mass(5);
    double mass7=mass(6);
    dxdt[0] = x[1]; // velocity in x1
    dxdt[1] = (newtonsg*(((mass2*(x[6]-x[0]))/(denom1))+(mass3*(x[12]-x[0])/(denom2))+(mass4*(x[18]-x[0])/denom4)+(mass5*(x[24]-x[0])/denom7)+(mass6*(x[30]-x[0])/denom11)+(mass7*(x[36]-x[0])/denom16)));
    dxdt[2] = x[3]; //velocity in y1
    dxdt[3] = (newtonsg*(mass2*(x[8]-x[2])/denom1+mass3*(x[14]-x[2])/denom2+(mass4*(x[20]-x[2])/denom4)+(mass5*(x[26]-x[2])/denom7)+(mass6*(x[32]-x[2])/denom11)+(mass7*(x[38]-x[2])/denom16)));
    dxdt[4] = x[5]; //velocity in z1
    dxdt[5] = (newtonsg*((mass2*(x[10]-x[4])/denom1)+(mass3*(x[16]-x[4])/denom2)+(mass4*(x[22]-x[4])/denom4)+(mass5*(x[28]-x[4])/denom7)+(mass6*(x[34]-x[4])/denom11)+(mass7*(x[40]-x[4])/denom16)));
    dxdt[6] = x[7]; //velocity in x2
    dxdt[7] =newtonsg*(mass1*(x[0]-x[6])/denom1+mass3*(x[12]-x[6])/denom3+(mass4*(x[18]-x[6])/denom5)+(mass5*(x[24]-x[6])/denom8)+(mass6*(x[30]-x[6])/denom12)+(mass7*(x[36]-x[6])/denom17));
    dxdt[8] = x[9]; //velocity in y2
    dxdt[9] =newtonsg*( (mass1*(x[2]-x[8])/denom1)+mass3*(x[14]-x[8])/denom3+(mass4*(x[20]-x[8])/denom5)+(mass5*(x[26]-x[8])/denom8)+(mass6*(x[32]-x[8])/denom12)+(mass7*(x[38]-x[8])/denom17));
    dxdt[10] = x[11]; //velocity in z2
    dxdt[11] =newtonsg*( (mass1*(x[4]-x[10])/denom1)+mass3*(x[16]-x[10])/denom3+(mass4*(x[22]-x[10])/denom5)+(mass5*(x[28]-x[10])/denom8)+(mass6*(x[34]-x[10])/denom12)+(mass7*(x[40]-x[10])/denom17));
    dxdt[12] = x[13]; //velocity in x3
    dxdt[13] =newtonsg*( (mass2*(x[6]-x[12])/denom3)+mass1*(x[0]-x[12])/denom2+(mass4*(x[18]-x[12])/denom6)+(mass5*(x[24]-x[12])/denom9)+(mass6*(x[30]-x[12])/denom13)+(mass7*(x[36]-x[12])/denom18));
    dxdt[14] = x[15]; //velocity in y3
    dxdt[15] =newtonsg*( (mass2*(x[8]-x[14])/denom3)+mass1*(x[2]-x[14])/denom2+(mass4*(x[20]-x[14])/denom6)+(mass5*(x[26]-x[14])/denom9)+(mass6*(x[32]-x[14])/denom13)+(mass7*(x[38]-x[14])/denom18));
    dxdt[16] = x[17]; //velocity in z3
    dxdt[17] =newtonsg*((mass2*(x[10]-x[16])/denom3)+mass1*(x[4]-x[16])/denom2+(mass4*(x[22]-x[16])/denom6)+(mass5*(x[28]-x[16])/denom9)+(mass6*(x[34]-x[16])/denom13)+(mass7*(x[40]-x[16])/denom18));
    dxdt[18]= x[19];//velocity x4
    dxdt[19]= newtonsg*((mass1*(x[0]-x[18])/denom4)+mass2*(x[6]-x[18])/denom5+(mass3*(x[12]-x[18])/denom6)+(mass5*(x[24]-x[18])/denom10)+(mass6*(x[30]-x[18])/denom14)+(mass7*(x[36]-x[18])/denom19));
    dxdt[20]= x[21];//velocity y4
    dxdt[21]= newtonsg*((mass1*(x[2]-x[20])/denom4)+mass2*(x[8]-x[20])/denom5+(mass3*(x[14]-x[20])/denom6)+(mass5*(x[26]-x[20])/denom10)+(mass6*(x[32]-x[20])/denom14)+(mass7*(x[38]-x[20])/denom19));
    dxdt[22]= x[23]; //velocity z4
    dxdt[23]= newtonsg*((mass1*(x[4]-x[22])/denom4)+mass2*(x[10]-x[22])/denom5+(mass3*(x[16]-x[22])/denom6)+(mass5*(x[28]-x[22])/denom10)+(mass6*(x[34]-x[22])/denom14)+(mass7*(x[40]-x[22])/denom19));
    dxdt[24]=x[25];
    dxdt[25]=newtonsg*((mass1*(x[0]-x[24])/denom7)+mass2*(x[6]-x[24])/denom8+(mass3*(x[12]-x[24])/denom9)+(mass4*(x[18]-x[24])/denom10)+(mass6*(x[30]-x[24])/denom15)+(mass7*(x[36]-x[24])/denom20));
    dxdt[26]=x[27];
    dxdt[27]=newtonsg*((mass1*(x[2]-x[26])/denom7)+mass2*(x[8]-x[26])/denom8+(mass3*(x[14]-x[26])/denom9)+(mass4*(x[20]-x[26])/denom10)+(mass6*(x[32]-x[26])/denom15)+(mass7*(x[38]-x[26])/denom20));
    dxdt[28]=x[29];
    dxdt[29]=newtonsg*((mass1*(x[4]-x[28])/denom7)+mass2*(x[10]-x[28])/denom8+(mass3*(x[16]-x[28])/denom9)+(mass4*(x[22]-x[28])/denom10)+(mass6*(x[34]-x[28])/denom15)+(mass7*(x[40]-x[28])/denom20));
    dxdt[30]=x[31];
    dxdt[31]=newtonsg*((mass1*(x[0]-x[30])/denom11)+mass2*(x[6]-x[30])/denom12+(mass3*(x[12]-x[30])/denom13)+(mass4*(x[18]-x[30])/denom14)+(mass5*(x[24]-x[30])/denom15)+(mass7*(x[36]-x[30])/denom21));
    dxdt[32]=x[33];
    dxdt[33]=newtonsg*((mass1*(x[2]-x[32])/denom11)+mass2*(x[8]-x[32])/denom12+(mass3*(x[14]-x[32])/denom13)+(mass4*(x[20]-x[32])/denom14)+(mass5*(x[26]-x[32])/denom15)+(mass7*(x[38]-x[32])/denom21));
    dxdt[34]=x[35];
    dxdt[35]=newtonsg*((mass1*(x[4]-x[34])/denom11)+mass2*(x[10]-x[34])/denom12+(mass3*(x[16]-x[34])/denom13)+(mass4*(x[22]-x[34])/denom14)+(mass5*(x[28]-x[34])/denom15)+(mass7*(x[40]-x[34])/denom21));
    dxdt[36]=x[37];
    dxdt[37]=newtonsg*((mass1*(x[0]-x[36])/denom16)+mass2*(x[6]-x[36])/denom17+(mass3*(x[12]-x[36])/denom18)+(mass4*(x[18]-x[36])/denom19)+(mass5*(x[24]-x[36])/denom20)+(mass6*(x[30]-x[36])/denom21));
    dxdt[38]=x[39];
    dxdt[39]=newtonsg*((mass1*(x[2]-x[38])/denom16)+mass2*(x[8]-x[38])/denom17+(mass3*(x[14]-x[38])/denom18)+(mass4*(x[20]-x[38])/denom19)+(mass5*(x[26]-x[38])/denom20)+(mass6*(x[32]-x[38])/denom21));
    dxdt[40]=x[41];
    dxdt[41]=newtonsg*((mass1*(x[4]-x[40])/denom16)+mass2*(x[10]-x[40])/denom17+(mass3*(x[16]-x[40])/denom18)+(mass4*(x[22]-x[40])/denom19)+(mass5*(x[28]-x[40])/denom20)+(mass6*(x[34]-x[40])/denom21));
}

void rhs_algorithm (const state_type &x, state_type &dxdt, const double t)
{       //this was the routine that was sipposed to solve all the derivatives

    double cup[number_bodies*6];
    for(int o=0;o<(number_bodies*3);o++)
        dxdt[(2*o)]=x[(o*2)+1];

    for(int i=0;i<number_bodies;i++)
    {
        int ib=i*6;
        for(int ic=1;ic<4;ic++)
        {
            cup[ib+(2*ic)-1]=0;
            for(int j=0;j<number_bodies;j++)
            {
                int jb=j*6;

                if(i!=j)
                {   double massb=mass(j);
                    double denom=pow(sqrt(get_denom(x,i,j)),3);
                    //cout<<denom<<endl;

                    cup[ib+(2*ic)-1]+=newtonsg*(massb*((x[ib+(2*ic)-2]-x[jb+(2*ic)-2])/denom));

                }
            }
            dxdt[ib+(2*ic)-1]=cup[ib+(2*ic)-1];
            std::cout<<dxdt[ib+(2*ic)-1]<<endl;

        }
    }
}

ofstream out;
void write_out( const state_type &x , const double t ) //outputs the data.
{
    out << t << '\t' << x[0] << '\t' << x[2]<<'\t'<< x[6] <<'\t'<< x[8] <<'\t'<<x[12]<<'\t'<<x[14]<<'\t'<< x[18]<<'\t'<<x[20]<<'\t'
    <<x[24]<<'\t'<<x[26]<<'\t'<<x[30]<<'\t'<<x[32]<<'\t'<<x[36]<<'\t'<<x[38]<<'\t'<<endl;// x[12]<<'\t'<<x[14]<<'\t'<< endl;
}

void compute_kinetic(const state_type &x , state_type_4  &kinetic_energy) //computes the kinetic energy of a body
{
    double mass_1=0.0;
    for(int i=0;i<number_bodies;i++)
    {
        int ib=6*i;
        mass_1 = mass(i);
        kinetic_energy[i]=(0.5)*mass_1*(pow(x[ib+1],2)+pow(x[ib+3],2)+pow(x[ib+5],2));
    }

}

void compute_potential(const state_type &x, state_type_3 &gravitational_energy)
{           //compute the potential energy of a body
    double mass_1 = mass(0);
    double mass_2 = mass(1);
    double mass_3 = mass(2);
    double mass_4 = mass(3);
    double mass_5 = mass(4);
    double mass_6 = mass(5);
    double mass_7 = mass(6);
    double distance_r_1=sqrt(get_denom(x,0,1)); // distance between sun and earth
    double distance_r_2=sqrt(get_denom(x,0,2));  // distance between sun and mars
    double distance_r_3=sqrt(get_denom(x,1,2));  // distance between earth and mars.
    double denom4=sqrt(get_denom(x,0,3)); // sun and jupiter
    double denom5=sqrt(get_denom(x,1,3)); // earth and jupiter
    double denom6=sqrt(get_denom(x,2,3));// mars and jupiter
    double denom7=sqrt(get_denom(x,0,4)); //sun and saturn
    double denom8=sqrt(get_denom(x,1,4)); //earth and saturn
    double denom9=sqrt(get_denom(x,2,4)); //mars and saturn
    double denom10=sqrt(get_denom(x,3,4)); //jupiter and saturn
    double denom11=sqrt(get_denom(x,0,5)); //sun and uranus
    double denom12=sqrt(get_denom(x,1,5)); //earth and uranus
    double denom13=sqrt(get_denom(x,2,5)); //mars and uranus
    double denom14=sqrt(get_denom(x,3,5)); //jupiter and uranus
    double denom15=sqrt(get_denom(x,4,5)); //saturn adn uranus
    double denom16=sqrt(get_denom(x,0,6)); //sun and neptune
    double denom17=sqrt(get_denom(x,1,6)); //earth and neptune
    double denom18=sqrt(get_denom(x,2,6)); //mars and neptune
    double denom19=sqrt(get_denom(x,3,6)); //jupiter and neptune
    double denom20=sqrt(get_denom(x,4,6)); //saturn and neptune
    double denom21=sqrt(get_denom(x,5,6));  //uranus and neptune

    gravitational_energy[0]=newtonsg*mass_1*((mass_2/distance_r_1)+(mass_3/distance_r_2)+(mass_4/denom4)+(mass_5/denom7)+(mass_6/denom11)+(mass_7/denom16));
    gravitational_energy[1]=newtonsg*mass_2*((mass_1/distance_r_1)+(mass_3/distance_r_3)+(mass_4/denom5)+(mass_5/denom8)+(mass_6/denom12)+(mass_7/denom17));
    gravitational_energy[2]=newtonsg*mass_3*((mass_1/distance_r_2)+(mass_2/distance_r_3)+(mass_4/denom6)+(mass_5/denom9)+(mass_6/denom13)+(mass_7/denom18));
    gravitational_energy[3]=newtonsg*mass_4*((mass_1/denom4)+(mass_2/denom5)+(mass_3/denom6)+(mass_5/denom10)+(mass_6/denom14)+(mass_7/denom19));
    gravitational_energy[4]=newtonsg*mass_5*((mass_1/denom7)+(mass_2/denom8)+(mass_3/denom9)+(mass_4/denom10)+(mass_6/denom15)+(mass_7/denom20));
    gravitational_energy[5]=newtonsg*mass_6*((mass_1/denom11)+(mass_2/denom12)+(mass_3/denom13)+(mass_4/denom14)+(mass_5/denom15)+(mass_7/denom21));
    gravitational_energy[6]=newtonsg*mass_7*((mass_1/denom16)+(mass_2/denom17)+(mass_3/denom18)+(mass_4/denom19)+(mass_5/denom20)+(mass_6/denom21));
}

void calculate_energy(const state_type &x, state_type_1 &error_array) //combines the potential adn kinetic enrgy
{
    state_type_4 gravitational_energy={0,0,0,0,0,0,0};
    error_array[0]=0.0;
    error_array[1]=0.0;
    error_array[2]=0.0;
    error_array[3]=0.0;
    error_array[4]=0.0;
    error_array[5]=0.0;
    error_array[6]=0.0;
    state_type_3 kinetic_energy={0,0,0,0,0,0,0};
    compute_kinetic(x, kinetic_energy);
    compute_potential(x, gravitational_energy);
    for(int i=0;i<number_bodies;i++)
            cout<<gravitational_energy[i]<<endl;
    for(int i=0;i<number_bodies;i++)
            cout<<kinetic_energy[i]<<endl;
    cout<<endl;


    for(int i=0;i<number_bodies;i++)
        error_array[i]=gravitational_energy[i]+kinetic_energy[i];
}

int main()
{
    bulirsch_stoer_dense_out< state_type > stepper( 1E-8 , 0.0 , 0.0 , 0.0 ); //different stepperd by BOOST
    bulirsch_stoer< state_type > stepper2( 1E-8 , 0.0 , 0.0 , 0.0 );

    state_type x = {{ 2.0 / sqrt(3.0),1,-10,1,-30,1,100,1,-50,1,50,1,40,1,20,1,-55,2,0,0,0,0,0,0,0,0,0,0}};

    double t = 0.00;
    //double t = 0.0;
    double dt = 1800;
    //double t_end = M_PI/2.0 - 0.1;
    double t_end = 193000000.0;

    //the following list is all of the initial conditions

    x[0] = -5000000.0; //# x1   //star
    x[1] = 0.0; //# vx1
    x[2] = 0.0 ;//# y1
    x[3] = -3.3 ; //vy1
    x[4] = -1.01 ;//# z1
    x[5] = 1.01 ;//# vz1

    x[6] = 1.361000e10*(0.980111);//# x2  //planet 1   kepler-11b
    x[7] =  9.6554e4*(0.19844);//# vx2
    x[8] = 1.361000e10*(-0.19844);//# y2
    x[9] =  9.6554e4*(0.980111);//# vy2
    x[10] = -6.0 ;//# z2
    x[11] = 1.000;//# vz2

    x[12] =  1.5993255e10;//# x3  //planet 2 kepler-11 c
    x[13] = 0.0;//# vx3
    x[14] =  0.0;//# y3
    x[15] =  8.9301e4;//# vy3
    x[16] = 0.000 ;//# z3
    x[17] = 0.000 ;//# vz3

    x[18]= (-0.95687)*2.315429e10; //x4  //planet 3 kepler-11 d
    x[19]= 7.4219e4*(0.290512);  //vx4
    x[20]= (-0.290512)*2.315429e10;  //y4
    x[21]= 7.4219e4*(-0.95687);  //vy4
    x[22]= 0.0; //z4
    x[23]= 0.0;  //vz4

    x[24]= 2.91187e10; //x5
    x[25]= (0.0025889)*6.6184e4; //vx5 planet 4 kepler-11 e
    x[26]= -(0.0025889)*2.91187e10; //y5
    x[27]= 6.6184e4;//vy5
    x[28]= 1.0; //z5
    x[29]= 1.0;  //vz5

    x[30]=3.7399000e10*(0.6440731); // planet 5 kepler-11 f
    x[31]=(-0.764963)*5.8252000e4;
    x[32]=3.7399000e10*(0.7649639);
    x[33]=5.8252000e4*(0.6440731);
    x[34]=1.0;
    x[35]=1.0;

    x[36]=6.97126000e10*(-0.060356); //planet 6 kepler-11 g
    x[37]=4.2824000e4*(0.998176);
    x[38]=6.97126000e10*(-0.998176);
    x[39]=-4.2824000e4*(0.0660356);
    x[40]=1.0;
    x[41]=1.0;

    state_type_1 error_array={{0,0,0,0,0,0,0}};     //this deal with the energy conservation
    state_type_2 error_start={{0,0,0,0,0,0,0}};
    calculate_energy(x,error_array);
    for(int y=0;y<number_bodies;y++)
        error_start[y]=error_array[y];
    double total_error_start=error_start[0]+error_start[1]+error_start[2]+error_start[3]+error_start[4]+error_start[5]+error_start[6];


    out.open( "bs3.dat" ); //fils the spreadsheet
    out.precision(16);
    integrate_adaptive( stepper , rhs , x , t , t_end , dt , write_out );
    out.close();



    calculate_energy(x,error_array); //final energy summation

    double total_error_finish=error_array[0]+error_array[1]+error_array[2]+error_array[3]+error_array[4]+error_array[5]+error_array[6];
    cout<<endl;
    std::cout<<(total_error_start/total_error_finish)<<endl;
    cout<<total_error_finish<<'\t'<<total_error_start<<endl;

    typedef runge_kutta_dopri5< state_type > dopri5_type;
    typedef controlled_runge_kutta< dopri5_type > controlled_dopri5_type;
    typedef dense_output_runge_kutta< controlled_dopri5_type > dense_output_dopri5_type;

    dense_output_dopri5_type dopri5 = make_dense_output( 1E-9 , 1E-9 , dopri5_type() );


}
