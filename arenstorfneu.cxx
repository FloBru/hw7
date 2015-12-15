#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

void func(double* k);
void errors(double* y, double* y5, double& error, double& maxerror);

int main() {
    
    double y[4]; //rk4 //y=(x,x`,y,y`) 
    double y5[4]; //rk5
    y[0] = 0.994;   //Startwerte rk4
    y[1] = 0.; 
    y[2] = 0.;
    y[3] = (-1)*2.00158510637908;
    y5[0] = 0.994;   //Startwerte rk5
    y5[1] = 0.; 
    y5[2] = 0.;
    y5[3] = (-1)*2.00158510637908;
    
    double time =0;
    int tend = 1; //Zeitendwert
    double dt = 0.1; //Schrittweite
    double error = 0.;
    double maxerror = 0.;
    double tol = pow(10., -5);
    
    
    double k1[4];  //RK4: 7k mit 4 Dimensionen
    double k2[4];
    double k3[4];
    double k4[4];
    double k5[4];
    double k6[4];
    double k7[4];
    
    double k51[4];  //RK5: 7k mit 4 Dimensionen
    double k52[4];
    double k53[4];
    double k54[4];
    double k55[4];
    double k56[4];
    double k57[4];
    
    ofstream out("arenstorf.txt");
    out << 0 << "\t" << y[0] << "\t" << y[1] << "\t" << y[2] << y[3] << endl;
    
    while(time <= (tend -dt)) {
        time += dt; 
        //k1
        for(int i = 0; i < 4; i++) {
            k1[i] = y[i];
            k51[i] = y5[i];
        }
        func(k1);
        func(k51);
    
        //k2
        for(int i = 0; i < 4; i++) {
            k2[i] = y[i] + (dt/5)*k1[i];
            k52[i] = y5[i] + (dt/5)*k51[i];
        }
        func(k2);
        func(k52);

        //k3
        for(int i = 0; i < 4; i++) {
            k3[i] = y[i] + dt*((3./40)*k1[i] + (9./40)*k2[i]);
            k53[i] = y5[i] + dt*((3./40)*k51[i] + (9./40)*k52[i]);
        }
        func(k3);
        func(k53);
        
        //k4
        for(int i = 0; i < 4; i++) {
            k4[i] = y[i] + dt*((44./45)*k1[i] + (-(56.)/15)*k2[i] + (32./9)*k3[i]);
            k54[i] = y5[i] + dt*((44./45)*k51[i] + (-(56.)/15)*k52[i] + (32./9)*k53[i]);
        }
        func(k4);
        func(k54);
        
        //k5
        for(int i = 0; i < 4; i++) {
            k5[i] = y[i] + dt*((19372./6561)*k1[i] + (-(25360.)/2187)*k2[i] + (64448./6561)*k3[i] + (-(212.)/729)*k4[i]);
            k55[i] = y5[i] + dt*((19372./6561)*k51[i] + (-(25360.)/2187)*k52[i] + (64448./6561)*k53[i] + (-(212.)/729)*k54[i]);
        }
        func(k5);
        func(k55);
        
        //k6
        for(int i = 0; i < 4; i++) {
            k6[i] = y[i] + dt*((9017./3168)*k1[i] + (-(355.)/33)*k2[i] + (46732./5247)*k3[i] + (49./176)*k4[i] + (-(5103.)/18656)*k5[i]);
            k56[i] = y5[i] + dt*((9017./3168)*k51[i] + (-(355.)/33)*k52[i] + (46732./5247)*k53[i] + (49./176)*k54[i] + (-(5103.)/18656)*k55[i]);
        }
        func(k6);
        func(k56);
        
        //k7
        for(int i = 0; i < 4; i++) {
            k7[i] = y[i] + dt*((35./384)*k1[i] + (500./1113)*k3[i] + (125./192)*k4[i] + (-(2187.)/6784)*k5[i] + (11./84)*k6[i]);
            k57[i] = y5[i] + dt*((35./384)*k51[i] + (500./1113)*k53[i] + (125./192)*k54[i] + (-(2187.)/6784)*k55[i] + (11./84)*k56[i]);
        }
        func(k7);
        func(k57);
        
        //yn+1 rk4
        for(int j = 0; j < 4; j++) {
            y[j] += dt*((35./384)*k1[j] + (500./1113)*k3[j] + (125./192)*k4[j] + (-2187/6784)*k5[j] + (11/84)*k6[j]);
        }
        
        //yn+1 rk5
        for(int j = 0; j < 4; j++) {
            y5[j] += dt*((5179./57600)*k51[j] + (7571./16695)*k53[j] + (393./640)*k54[j] + (-92097/339200)*k55[j] + (187/2100)*k56[j] + (1/40)*k57[j]);
        }
        
        errors(y, y5, error, maxerror);
        //wenn der Fehler zu gross ist, dann wird eine neue schrittweite gewaehlt und mit rk4 weitergerechnet, q=0,5, p=4
        if (maxerror > tol) {
            dt *= 0.5*pow(tol/maxerror, 1./(4+1));
            for (int i = 0; i < 4; i++) {
                y5[i] = y[i];
            }
        }
        out << time << "\t" << y[0] << y [1] << "\t" << y[2] << y [3] << "\t" << endl;
            
    }
    out.close();
    return 0;
}
//arenstorf-orbits
void func(double* k){
    double yt[4];
    const double mass = 0.012277471;
    double r = sqrt((yt[0]+yt[2])*(yt[0]+yt[2]) + yt[2]*yt[2]);
    double s = sqrt((yt[0]-1-mass)*(yt[0]-1-mass) + yt[2]*yt[2]);
     
    for(int i = 0; i < 4; i++) {
        yt[i] = k[i];
    }
    k[0] = yt[1];
    k[1] = yt[0] + 2*yt[2] - ((1-mass)*(yt[0]+mass))/(r*r*r) - (mass*(yt[0]-1+mass))/(s*s*s);
    k[2] = yt[3];
    k[3] = yt[2] - 2*yt[1] - ((1-mass)*yt[2])/(r*r*r) - (mass*yt[2]/(s*s*s));
}
//Fehler
void errors(double* y, double* y5, double& error, double& maxerror) {
    maxerror = 0.;
    for(int i; i < 4; i++) {
        error = abs(y[i] - y5[i]);
        if(error > maxerror) {
            maxerror = error;
        }
    }
}
        
        
        
        
        
        
        
        