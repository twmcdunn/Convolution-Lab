import java.util.ArrayList;
import java.util.Arrays;
import java.awt.Graphics;
import java.awt.Color;
import javax.swing.*;
/**
3d rep of a sound in: time, freq, amp
 */
public class Sound3D
{
    public double fund;
    public ArrayList<ArrayList<double[]>> mySpects;
    public static int numOfSamples = 30 * 20;
    public static int windowSize = (int)Math.pow(2,11);
    public ArrayList<double[]> aveSpect;
    public ArrayList<Integer> peakFreqs;
    public ArrayList<Envelope> peakEnvs;
    public Sound3D(String fileName)
    {
        mySpects = new ArrayList<ArrayList<double[]>>();
        loadSpects(fileName);
        //display();
        System.out.println("Spects Loaded");
    }

    public Sound3D(Sound3D s1,Sound3D s2){
        mySpects = new ArrayList<ArrayList<double[]>>();
        convolution(s1,s2);
    }

    private void convolution(Sound3D s1,Sound3D s2){
        if(s1.fund > s2.fund){
            Sound3D s3 = s1;
            s1 = s2;
            s2 = s3;
        }
        fund = s1.fund;
        for(int n = 0; n < numOfSamples; n++){
            ArrayList<double[]> spect1 = s1.mySpects.get(n);
            ArrayList<double[]> spect2 = s2.mySpects.get(n);
            ArrayList<double[]> spectProduct = new ArrayList<double[]>();
            for(int i = 0; i < spect1.size(); i++){
                double exInd = i * s2.fund / s1.fund;
                int index = (int)exInd;
                double fract = exInd - index;
                double amp2 = 0;

                if(index < spect2.size() - 2){
                    amp2 = spect2.get(index)[0] * (1-fract) + spect2.get(index)[1] * fract;
                }
                spectProduct.add(new double[]{spect1.get(i)[0],spect1.get(i)[1]*amp2});
            }
            mySpects.add(spectProduct);
        }
        loadAveSpect();
    }

    private void loadSpects(String sample){

        float[] read = ReadSound.readSound(sample);
        double[] sound = new double[read.length];
        for(int i = 0; i< sound.length; i++){
            sound[i] = read[i];
        }

        for(int a = 0; a < numOfSamples; a++){
            int stWind = (int)(a * (sound.length - windowSize) / (double) numOfSamples);

            double[] samp = Arrays.copyOfRange(sound, stWind, stWind + windowSize);
            for(int i = 0; i < samp.length; i++){
                double windowFunction = (2 - (Math.cos(Math.PI * 2 *  i / (double)(windowSize - 1)) + 1)) / 2.0;
                samp[i] *= windowFunction;
            }

            ArrayList<double[]> spec = FFT.getSpectrumWPhase(samp, false);
            mySpects.add(spec);
        }
        // display();
        double maxVol = Double.MIN_VALUE;
        for(ArrayList<double[]> spect: mySpects)
            for(double[] bin: spect)
                maxVol = Math.max(maxVol, bin[1]);
        for(ArrayList<double[]> spect: mySpects)
            for(double[] bin: spect)
                bin[1] /= maxVol;
        loadAveSpect();
        loadPeakEnvs();

        maxVol = 0;//now used for ave of all samples
        for(int i = 0; i < aveSpect.size(); i++){
            if(aveSpect.get(i)[1] > maxVol){
                maxVol = aveSpect.get(i)[1];
                fund = aveSpect.get(i)[0] * WaveWriter.SAMPLE_RATE;
            }
        }
    }

    public void loadAveSpect(){
        aveSpect = new ArrayList<double[]>();
        for(int i = 0; i < mySpects.get(0).size(); i++)
            aveSpect.add(new double[2]);
        for(ArrayList<double[]> spect: mySpects)
            for(int i = 0; i < spect.size(); i++){
                aveSpect.get(i)[0] = spect.get(i)[0];
                aveSpect.get(i)[1] += spect.get(i)[1] / (double) numOfSamples;
            }
        double maxVol = 0;
        for(double[] bin: aveSpect){
            maxVol = Math.max(maxVol, bin[1]);
        }
        for(double[] bin: aveSpect){
            bin[1] /= maxVol;
        }
        int searchWindow = 3;
        peakFreqs = new ArrayList<Integer>();
        for(int i = 1; i < aveSpect.size() - 1; i++){
            if(aveSpect.get(i-1)[1] < 0.001)//0.1
                continue;
            int n = 0;
            do{
                n++;
            }while(n <= searchWindow
            && aveSpect.get(i)[1] > aveSpect.get(i - n)[1]
            && aveSpect.get(i)[1] > aveSpect.get(i + n)[1]);

            if(n > searchWindow){
                peakFreqs.add(i);
            }
            /*
            if(aveSpect.get(i-1)[1] < aveSpect.get(i)[1] 
            && aveSpect.get(i+1)[1] < aveSpect.get(i)[1]){
            peakFreqs.add(i);
            }
             */
        }
        System.out.println("DETECTED PEAKS: ");
        for(int i = 0; i < peakFreqs.size(); i++){
            System.out.println(mySpects.get(0).get(peakFreqs.get(i))[0] * WaveWriter.SAMPLE_RATE);
        }
    }

    public void loadPeakEnvs(){
        peakEnvs = new ArrayList<Envelope>();
        for(int n: peakFreqs){
            double[] vals = new double[numOfSamples];
            double[] times = new double[numOfSamples];
            for(int i = 0; i < numOfSamples; i++){
                times[i] = i / (double)(numOfSamples - 1);
                vals[i] = mySpects.get(i).get(n)[1];
            }
            peakEnvs.add(new Envelope(times, vals));
        }
    }

    double minFreq = Double.MAX_VALUE;
    double maxFreq = Double.MIN_VALUE;
    private void display(){
        for(ArrayList<double[]> spect: mySpects)
            for(int n: peakFreqs){
                double[] bin = spect.get(n);
                minFreq = Math.min(bin[0], minFreq);
                maxFreq = Math.max(bin[0], maxFreq);
            }
        JPanel jp = new JPanel(){
                // double[] myNormVals;
                @Override
                public void paint(Graphics g){
                    g.setColor(Color.WHITE);
                    g.fillRect(0,0,500,500);

                    g.setColor(Color.BLACK);
                    for(int i = 0; i < mySpects.size(); i++){
                        for(int n: peakFreqs){
                            double d = (mySpects.get(i).get(n)[0] - minFreq) / (maxFreq - minFreq);
                            d *= 500;
                            int x = (int)(500 *i / (double)mySpects.size());
                            double amp = mySpects.get(i).get(n)[1];
                            g.setColor(new Color((int)(255 * amp),(int)(255 * amp),(int)(255 * amp)));
                            //g.setColor(Color.BLACK);
                            g.fillOval( x, (int)(500 -  d), 1, 1);
                        }
                    }

                }
            };
        JFrame jf = new JFrame("DELTA VALUES");
        jf.setBounds(0,0,500,500);
        jf.add(jp);
        jf.setVisible(true);
        jp.repaint();

    }

    int grainSize = 5000;
    public double[] tone(double freq, double vol, double dur){
        double[] tone = new double[(int)(dur*WaveWriter.SAMPLE_RATE)];
        for(int a = 0; a < 100000; a++){
            int n = (int)(mySpects.get(0).size() * Math.random());
            double f = mySpects.get(0).get(n)[0] * WaveWriter.SAMPLE_RATE * freq / fund;
            //System.out.println("BF: " + f);
            double t = Math.random();
            double exInd = (numOfSamples - 1) * t;
            int index = (int)Math.rint(exInd);
            //double fract = exInd - index;
            double amp = mySpects.get(index).get(n)[1]; //= mySpects.get(index).get(n)[1] * (1-fract) + mySpects.get(index + 1).get(n)[1] * fract;
            for(int i = 0; i < grainSize; i++){
                double windowFunction = (2 - (Math.cos(Math.PI * 2 *  i / (double)(grainSize - 1)) + 1)) / 2.0;
                tone[i + (int)(t * (tone.length - grainSize))] += windowFunction * amp * Math.sin(Math.PI * 2 * f * i / (double) WaveWriter.SAMPLE_RATE);
            }
        }
        return tone;
    }

    public double[] tone1(double freq, double vol, double dur){
        double[] tone = new double[(int)(dur*WaveWriter.SAMPLE_RATE)];
        for(int n : peakFreqs){
            double f = mySpects.get(0).get(n)[0] * WaveWriter.SAMPLE_RATE * freq / fund;
            for(int i = 0; i < tone.length; i++){
                double t = i / (double)tone.length;
                double exInd = (numOfSamples - 1) * t;
                int index = (int)exInd;
                double fract = exInd - index;
                double amp = mySpects.get(index).get(n)[1] * (1-fract) + mySpects.get(index + 1).get(n)[1] * fract;;

                tone[i] +=  amp * Math.sin(Math.PI * 2 * f * i / (double) WaveWriter.SAMPLE_RATE);
            }
        }
        return tone;
    }

    public double[] tone2(double freq, double vol, double dur){
        grainLength = WaveWriter.SAMPLE_RATE / 2;
        double[] tone = new double[(int)(dur*WaveWriter.SAMPLE_RATE)];
        //int grainLength = 2 * (int)(tone.length / (double)(numOfSamples + 2));
        for(int samp = 0; samp < numOfSamples; samp++){
            for(int n : peakFreqs){
                double f = mySpects.get(0).get(n)[0] * WaveWriter.SAMPLE_RATE * freq / fund;
                double amp = mySpects.get(samp).get(n)[1];
                int offset = (int)(Math.random() * tone.length / (double) numOfSamples);
                for(int i = 0; i < grainLength 
                && i + (int)(tone.length * samp /(double) numOfSamples) + offset < tone.length; i++){
                    double windowFunction = (2 - (Math.cos(Math.PI * 2 *  i / (double)(grainLength - 1)) + 1)) / 2.0;
                    tone[i + (int)(tone.length * samp /(double) numOfSamples) + offset] += windowFunction * amp * Math.sin(Math.PI * 2 * f * i / (double) WaveWriter.SAMPLE_RATE);
                }
            }
        }
        return tone;
    }

    int grainLength = WaveWriter.SAMPLE_RATE / 10;
    public double[] tone3(double freq, double vol, double dur){
        double[] tone = new double[(int)(dur * WaveWriter.SAMPLE_RATE)];
        for(int a = 0; a < dur * 10000; a++){
            double t = Math.random();
            double exactInd = (t * numOfSamples) - 0.5;
            int index = 0;
            if(exactInd > numOfSamples - 1){
                index = numOfSamples - 1;
            }
            else if(exactInd < 0){
                index = 0;
            }
            else{
                index = (int)exactInd;
                double fract = exactInd - index;
                if(Math.random() < fract)
                    index++;
            }
            //System.out.println("SAMPLE INDEX:" + index);
            int freqInd = peakFreqs.get((int)(Math.random() * peakFreqs.size()));
            //freqInd = peakFreqs.get(0);
            double f = mySpects.get(index).get(freqInd)[0] * WaveWriter.SAMPLE_RATE * freq / fund;
            double amp = mySpects.get(index).get(freqInd)[1] * vol;
            //System.out.println("F: " + f);
            for(int i = 0; i < grainLength && i + (int)(t * tone.length) < tone.length; i++){
                double windowFunction = (2 - (Math.cos(Math.PI * 2 *  i / (double)(grainLength - 1)) + 1)) / 2.0;
                tone[i + (int)(t * tone.length)] += amp * windowFunction * Math.sin(2 * Math.PI * f * i / (double)WaveWriter.SAMPLE_RATE);
            }
            //break;
        }
        return tone;
    }

    public double[] tone4(double freq, double vol, double dur){
        double[] tone = new double[(int)(dur * WaveWriter.SAMPLE_RATE)];
        for(int a = 0; a < peakFreqs.size(); a++){
            int n = peakFreqs.get(a);
            double[] partial = new double[tone.length];
            for(int i = 0; i < partial.length; i++){
                partial[i] = Math.random() * 2 - 1;
            }
            for(int i = 0; i < 3; i++)
                partial = BPF.BPF(partial, WaveWriter.SAMPLE_RATE, mySpects.get(0).get(n)[0] * freq / fund, 30);
            for(int i = 0; i < partial.length; i++)
                partial[i] *= peakEnvs.get(a).getValue(i / (double)partial.length);
        }
        return tone;
    }

    public double[] tone5(double freq, double vol){
        grainSize = 2048;
        double[] tone = new double[mySpects.size() * grainSize / 2 + grainSize];
        for(int i = 0; i < mySpects.size(); i++){
            for(int n = 0; n < mySpects.get(i).size(); n++){
                double phase = Math.random() * Math.PI * 2;
                double f = mySpects.get(i).get(n)[0] * WaveWriter.SAMPLE_RATE * freq / fund;
                for(int j = 0; j < grainSize; j++){
                    double windowFunction = (2 - (Math.cos(Math.PI * 2 *  j / (double)(grainSize - 1)) + 1)) / 2.0;
                    tone[i*grainSize / 2 + j] += Math.sin(phase + Math.PI * 2 * f * j / (double)WaveWriter.SAMPLE_RATE);
                }
            }
        }
        return tone;
    }

    public static void test(){
        Sound3D a = new Sound3D("ss.wav");
        double[] tone = a.tone5(170, 1);//a.twoDInterp(170, 1, 3);//a.tone1(170 * 3, 1);

        WaveWriter ww = new WaveWriter("sound3dTest");

        for(int i = 0; i < tone.length; i++){
            ww.df[0][i] += tone[i];
            ww.df[1][i] += tone[i];
        }
        //System.out.println(min + " " + max);

        ww.render();
    }
}
