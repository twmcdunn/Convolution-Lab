import java.util.ArrayList;
import java.util.Arrays;
import javax.swing.*;
import java.awt.Graphics;
import java.awt.Color;
/**
 * Write a description of class Analysis here.
 *
 * @author (your name)
 * @version (a version number or a date)
 */
public class Analysis
{
    public String sample;
    public double fund;
    public int windowSize, fRes = 168;
    public ArrayList<ArrayList<double[]>> mySpects;
    public Analysis(String sampName){
        sample = sampName;
        fund = 170;
        windowSize = (int)Math.pow(2,11);//16);//11);
        mySpects = new ArrayList<ArrayList<double[]>>();
        loadSpect2();
    }

    public double[] tone(double freq, double vol){
        double[] tone = new double[(1 + mySpects.size()) * windowSize / 2];
        int extraBuff = (int)(WaveWriter.SAMPLE_RATE * 3);
        for(int i = 0; i < mySpects.size(); i++){
            for(int a = 0; a < mySpects.get(i).size(); a++){
                //a = 1;
                double[] bin = mySpects.get(i).get(a);
                double[] grain = new double[windowSize + extraBuff];
                for(int n = 0; n < grain.length; n++)
                    grain[n] = (Math.random() * 2 - 1);//*vol
                for(int n = 0; n < 5; n++)
                    grain = BPF.BPF(grain, WaveWriter.SAMPLE_RATE, bin[0] * freq, 1);
                grain = normalize(grain);
                for(int n = 0; n < windowSize; n++){
                    double windowFunction = (2 - (Math.cos(Math.PI * 2 *  n / (double)(windowSize - 1)) + 1)) / 2.0;
                    grain[n + extraBuff] *= windowFunction * bin[1];
                }
                for(int n = 0; n < windowSize; n++){
                    tone[n + i * windowSize / 2] += grain[n + extraBuff];
                }
                //break;
            }
            //break;
        }
        return tone;
    }

    public double[] tone1(double freq, double vol){
        double[] tone = new double[(1 + mySpects.size()) * windowSize / 2 + windowSize];
        for(int i = 0; i < mySpects.size(); i++){
            for(int a = 0; a < mySpects.get(i).size(); a++){
                // a = 1;
                double[] bin = mySpects.get(i).get(a);
                double[] grain = new double[windowSize];
                for(int b = 0; b < 1; b++){
                    double phase = Math.PI * 2;
                    double f = freq * bin[0];
                    for(int n = 0; n < grain.length; n++){
                        double windowFunction = (2 - (Math.cos(Math.PI * 2 *  n / (double)(windowSize - 1)) + 1)) / 2.0;
                        grain[n] = windowFunction * vol * bin[1] * Math.sin(phase + f * Math.PI * 2 * n / (double)(WaveWriter.SAMPLE_RATE));

                    }
                    int offSet = (int)(0.5 * windowSize * Math.random());

                    //offSet = 0;
                    for(int n = 0; n < windowSize; n++){
                        tone[offSet + n + i * windowSize / 2] += grain[n];
                    }
                }
                //break;
            }
            //break;
        }
        return tone;
    }

    public double[] twoDInterp(double freq, double vol, double dur){
        double[] tone = new double[(int)(dur * WaveWriter.SAMPLE_RATE)];
        for(int i = 0; i < 100000; i++){
            double t = Math.random();
            double[] grain = new double[WaveWriter.SAMPLE_RATE / 10];
            double[] freqAmp = chooseFreqAmp(t);
            double f = freq * freqAmp[0] / fund;
            double a = vol * freqAmp[1] ;
            for(int n = 0; n < grain.length; n++){
                double windowFunction = (2 - (Math.cos(Math.PI * 2 *  n / (double)(grain.length - 1)) + 1)) / 2.0;
                tone[(int)(t * (tone.length - grain.length)) + n] += a * windowFunction * Math.sin(f * 2 * Math.PI * n / (double)WaveWriter.SAMPLE_RATE);
            }

        }
        return tone;
    }

    private double[] chooseFreqAmp(double time){
        double exInd = ((mySpects.size() - 1) * time);
        int index = (int)exInd;
        double fract = exInd - index;
        ArrayList<double[]> firstSpect = mySpects.get(index);
        ArrayList<double[]> secondSpect = mySpects.get(index + 1);

        double[] vArr = new double[fRes];
        double totAmp = 0;
        for(int i = 0; i < fRes; i++){
            try{
                double freq = 20 * Math.pow(2,7 * i/(double)fRes);

                int fInd = 0;
                for(fInd = 0; fInd < firstSpect.size() && firstSpect.get(fInd)[0] > freq; fInd++){}
                double m = (firstSpect.get(fInd + 1)[1] - firstSpect.get(fInd)[1]) / (firstSpect.get(fInd + 1)[0] - firstSpect.get(fInd)[0]);
                double b = firstSpect.get(fInd)[1] - m * firstSpect.get(fInd)[0];
                double amp1 = freq * m + b;

                fInd = 0;
                for(fInd = 0; fInd < secondSpect.size() && secondSpect.get(fInd)[0] > freq; fInd++){}
                m = (secondSpect.get(fInd + 1)[1] - secondSpect.get(fInd)[1]) / (secondSpect.get(fInd + 1)[0] - secondSpect.get(fInd)[0]);
                b = secondSpect.get(fInd)[1] - m * secondSpect.get(fInd)[0];
                double amp2 = freq * m + b;

                double amp = amp2 * fract + amp1 * (1-fract);

                vArr[i] = amp;

                totAmp += amp;
            }
            catch(Exception e){}
        }
        double rand = Math.random();
        double totElapsed = 0;
        double[] freqAmp = null;
        for(int i = 0; i < vArr.length; i++){
            totElapsed += vArr[i] / totAmp;
            if(totElapsed > rand){
                double freq = 20 * Math.pow(2,7 * i/(double)fRes);
                double amp = vArr[i];
                freqAmp = new double[]{freq, amp};
                break;
            }
        }

        /*
        for(double v: vArr)
        System.out.println(v);
         */
        return freqAmp;
    }

    private static double[] normalize(double[] sig){
        double max = Double.MIN_VALUE;
        for(int i = 0; i < sig.length; i++){
            max = Math.max(max, Math.abs(sig[i]));
        }

        for(int i = 0; i < sig.length; i++){
            sig[i] /= max;
        }
        return sig;
    }

    public void loadSpect(){

        float[] read = ReadSound.readSound(sample);
        double[] sound = new double[read.length];
        for(int i = 0; i< sound.length; i++){
            sound[i] = read[i];
        }

        //14
        for(int stWind = 0; stWind < sound.length - windowSize; stWind += windowSize / 2){
            System.out.println(stWind / (double)WaveWriter.SAMPLE_RATE);
            double[] samp = Arrays.copyOfRange(sound, stWind, stWind + windowSize);
            for(int i = 0; i < samp.length; i++){
                double windowFunction = (2 - (Math.cos(Math.PI * 2 *  i / (double)(windowSize - 1)) + 1)) / 2.0;
                samp[i] *= windowFunction;
            }

            ArrayList<double[]> spec = FFT.getSpectrumWPhase(samp, false);
            ArrayList<double[]> spect1 = new ArrayList<double[]>();

            int searchWindow = 1;
            for(int i = searchWindow; i < spec.size() - searchWindow; i++){
                if(spec.get(i)[1] < 0.1)
                    continue;
                int n = 0;
                do{
                    n++;
                }while(n <= searchWindow
                && spec.get(i)[1] > spec.get(i - n)[1]
                && spec.get(i)[1] > spec.get(i + n)[1]);

                if(n > searchWindow){
                    spect1.add(new double[]{spec.get(i)[0] * WaveWriter.SAMPLE_RATE, spec.get(i)[1]});
                }
            }
            if(spect1.size() == 0){
                mySpects.add(spect1);
                continue;
            }
            if(fund == -1){
                fund = spect1.get(0)[0];
                fund = 170;
            }

            for(int i = 0; i < spect1.size(); i++){
                double freq = spect1.get(i)[0] / fund;
                spect1.get(i)[0] = freq;
            }
            mySpects.add(spect1);
        }

        double maxVol = Double.MIN_VALUE;
        for(ArrayList<double[]> spect: mySpects)
            for(double[] bin: spect)
                maxVol = Math.max(maxVol, bin[1]);
        for(ArrayList<double[]> spect: mySpects)
            for(double[] bin: spect)
                bin[1] /= maxVol;

        display();
    }

    public void loadSpect1(){
        float[] read = ReadSound.readSound(sample);
        double[] sound = new double[read.length];
        for(int i = 0; i< sound.length; i++){
            sound[i] = read[i];
        }

        //14
        for(int stWind = 0; stWind < sound.length - windowSize; stWind += windowSize / 2){
            System.out.println(stWind / (double)WaveWriter.SAMPLE_RATE);
            double[] samp = Arrays.copyOfRange(sound, stWind, stWind + windowSize);
            for(int i = 0; i < samp.length; i++){
                double windowFunction = (2 - (Math.cos(Math.PI * 2 *  i / (double)(windowSize - 1)) + 1)) / 2.0;
                samp[i] *= windowFunction;
            }

            ArrayList<double[]> spec = FFT.getSpectrumWPhase(samp, false);
            ArrayList<double[]> spect1 = new ArrayList<double[]>();

            int searchWindow = 1;
            for(int i = 0; i < spec.size(); i++){
                spect1.add(new double[]{spec.get(i)[0] * WaveWriter.SAMPLE_RATE, spec.get(i)[1]});
            }
            if(spect1.size() == 0){
                mySpects.add(spect1);
                continue;
            }
            if(fund == -1){
                fund = spect1.get(0)[0];
                fund = 170;
            }

            for(int i = 0; i < spect1.size(); i++){
                double freq = spect1.get(i)[0] / fund;
                spect1.get(i)[0] = freq;
            }
            mySpects.add(spect1);
        }

        double maxVol = Double.MIN_VALUE;
        for(ArrayList<double[]> spect: mySpects)
            for(double[] bin: spect)
                maxVol = Math.max(maxVol, bin[1]);
        for(ArrayList<double[]> spect: mySpects)
            for(double[] bin: spect)
                bin[1] /= maxVol;

        display();
    }

    public void loadSpect2(){

        float[] read = ReadSound.readSound(sample);
        double[] sound = new double[read.length];
        for(int i = 0; i< sound.length; i++){
            sound[i] = read[i];
        }

        //14
        for(int stWind = 0; stWind < sound.length - windowSize; stWind += windowSize / 2){
            System.out.println(stWind / (double)WaveWriter.SAMPLE_RATE);
            double[] samp = Arrays.copyOfRange(sound, stWind, stWind + windowSize);
            for(int i = 0; i < samp.length; i++){
                double windowFunction = (2 - (Math.cos(Math.PI * 2 *  i / (double)(windowSize - 1)) + 1)) / 2.0;
                samp[i] *= windowFunction;
            }

            ArrayList<double[]> spec = FFT.getSpectrumWPhase(samp, false);
            ArrayList<double[]> spect1 = new ArrayList<double[]>();

            int searchWindow = 1;
            for(int i = searchWindow; i < spec.size() - searchWindow; i++){
                if(spec.get(i)[1] < 0.1)
                    continue;
                int n = 0;
                do{
                    n++;
                }while(n <= searchWindow
                && spec.get(i)[1] > spec.get(i - n)[1]
                && spec.get(i)[1] > spec.get(i + n)[1]);

                if(n > searchWindow){
                    spect1.add(new double[]{spec.get(i)[0] * WaveWriter.SAMPLE_RATE, spec.get(i)[1]});
                }
            }
            if(spect1.size() == 0){
                mySpects.add(spect1);
                continue;
            }
            if(fund == -1){
                fund = spect1.get(0)[0];
                fund = 170;
            }

            for(int i = 0; i < spect1.size(); i++){
                double freq = spect1.get(i)[0] / fund;
                spect1.get(i)[0] = freq;
            }

            mySpects.add(spect1);
        }

        double maxVol = Double.MIN_VALUE;
        for(ArrayList<double[]> spect: mySpects)
            for(double[] bin: spect)
                maxVol = Math.max(maxVol, bin[1]);
        for(ArrayList<double[]> spect: mySpects)
            for(double[] bin: spect)
                bin[1] /= maxVol;

        display();
    }

    double minFreq = Double.MAX_VALUE;
    double maxFreq = Double.MIN_VALUE;
    public void display(){
        for(ArrayList<double[]> spect: mySpects)
            for(double[] bin: spect){
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
                        for(int n = 0; n < mySpects.get(i).size(); n++){
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

    public static void test(){
        Analysis a = new Analysis("ss.wav");
        double[] tone = a.tone1(170 * 3, 1);//a.twoDInterp(170, 1, 3);//a.tone1(170 * 3, 1);

        WaveWriter ww = new WaveWriter("resynthTest");

        for(int i = 0; i < tone.length; i++){
            ww.df[0][i] += tone[i];
            ww.df[1][i] += tone[i];
        }
        //System.out.println(min + " " + max);

        ww.render();
    }
}
