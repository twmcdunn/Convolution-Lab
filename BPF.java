import java.util.Arrays;
/**
 * Write a description of class BPF here.
 *
 * @author (your name)
 * @version (a version number or a date)
 */
public class BPF
{
    public static double[] BPF(double[] x, double sampleFreq, double centerFreq, double dur)
    {
        double[] y = new double[x.length];
        double q = (dur * centerFreq);// / A;
        double w0 = 2 * Math.PI * centerFreq / sampleFreq;
        double alfa = Math.sin(w0) / (2 * q);//Math.sin(w0) * Math.sinh((Math.log(2) / 2.0) * bw * (w0 / (Math.sin(20))));
        double b0 = q * alfa;//1 + A * alfa;
        double b1 = 0;//-2 * Math.cos(w0);
        double b2 = - q * alfa;;//1 - A * alfa;
        double a0 = 1 + alfa;
        double a1 = -2 * Math.cos(w0);
        double a2 = 1 - alfa;
        for(int n = 2; n < y.length; n++){
            y[n] = (b0/a0) * x[n] + (b1/a0) * x[n - 1] + (b2/a0) * x[n - 2] 
            - (a1/a0) * y[n - 1] - (a2/a0) * y[n - 2];
        }
        double onset = WaveWriter.SAMPLE_RATE/40.0;
        for(int i = 0; i < onset; i++){
            y[y.length - 1 - i] *= i / onset;
            y[i] *= i / onset;
        }
        return y;
    }


    public static void test1(){
        double[] sound = new double[WaveWriter.SAMPLE_RATE * 10];
        for(int i = 0; i < sound.length; i++){
            sound[i] = (Math.random() * 2 - 1)*0.01;
        }
        double[] noise = Arrays.copyOf(sound, sound.length);
        for(int i = 0; i < 3; i++)
            sound = BPF(Arrays.copyOf(sound, 60* WaveWriter.SAMPLE_RATE), WaveWriter.SAMPLE_RATE, 170, 1);
        for(int i = 0; i < noise.length; i++){
            sound[i] -= noise[i];
        }
        WaveWriter ww = new WaveWriter("test");

        for(int i = 0; i < sound.length; i++){
            ww.df[0][i] += sound[i];
            ww.df[1][i] += sound[i];
        }

        ww.render();

    }

}
