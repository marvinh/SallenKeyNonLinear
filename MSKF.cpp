//adapted from https://cytomic.com/files/dsp/SkfLinearTrapOptimised2.pdf 
//Convert to Linear Estimate, Tangential, then solve with Newton

    class MSKF {
    public:
        double sampleRate;
        MSKF() {
            
        }
        
        void setSampleRate(double newSampleRate) {
            sampleRate = newSampleRate;
            //onp.setParams(240, sampleRate); update one pole at 240hz  4.17ms
        }
        
        //onePoleLowpass onp; one pole lowpass smoother 
        
        
        double ic1eq = 0.0;
        double ic2eq = 0.0;
        double ic3eq = 0.0;
        
        void filter(double cutoff, double resonance, float *in, float * out, int blockSize) {
            //WARNGING: do not pass nyquist unless you oversample
            //CUTOFF = HZ/(SAMPLE_RATE*OVERSAMPLING_FACTOR)
            //RESONANCE = 0.0 - 1.0
            //std max will prevent self oscillation
            
            
            double g = tan(M_PI * cutoff);

            double k = (2+M_E)*resonance;
            double k2 = 2*resonance;
            
            //smooth cutoff because analog circuits are not discrete
            
            
            //g = onp.tickLinear(g);
            
            for (int i = 0; i < blockSize; i++) {
                
                //g = onp.tickLinear(g);
                
                double tg = 2*g;
                
                double g2 = g*g;
                
                double gd = 1.0/(-1 - tg - g2 + g*k2);
                
                double v0 = in[i];
                
                double error = 1;
                
                double tol = 1e-6;
                
                double v1 = 0.0;
                
                double v2 =  0.0;
                
                double gk = 1.0;
                               
                double g2v0 = g2*v0;
                double gv0 = g*v0;
                
                //LINEAR VOUT ESTIMATE
                v2 = -((g*ic1eq + ic1eq + g*ic2eq + g2v0)*gd);
                v1 = -(ic1eq - g*ic1eq - ic2eq*k2 - gv0 - g2v0)*gd;
                
                //Tangential Estimate
                //Use linear resonance to set tangential sounds better because we dont have ic3eq yet
                
                double base = k2*v2;
                float tBase = tanh(base);
                
                double a = sech2_with_tanh( tBase );
                double b = tBase - base;
                
                v2 = -(g*ic1eq + ic1eq + g*ic2eq + g2v0)/(-1.0 - tg - g2 + g*(a+b));
                
                
                double ggk = 1.0;
                double gic1eq = g*ic1eq;
                double gic2eq = g*ic2eq;
                double next = 0.0;
                for(int j = 0 ; j < 50; j++) {
                
                    double x = k * v2;
                    double fk = tanh(x);
                    gk = 1 - (fk*fk);
                    ggk = g*gk;
                    ic3eq = fk - gk*x;
                    next = -((gic1eq + ic2eq + gic2eq + g*ic3eq + g2v0)/(-1-tg-g2+ggk));
                    error = std::abs(next-v2);
                    if(error < tol){
                        v2 = next;
                        break;
                    }
                    v2 = next;
                }
                
                
                double omg = (-1-g);
                
                v1 = -((ic2eq * gk - omg * (ic1eq + ic3eq + gv0))/(-(omg*omg) + ggk));
                
                ic1eq = 2 * (v1 - (gk*v2 + ic3eq)) - ic1eq;
                
                ic2eq = 2 * v2 - ic2eq;
                
                out[i] = v2;
            }
            
        }
        
    };
