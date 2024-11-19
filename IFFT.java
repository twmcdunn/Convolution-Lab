import java.util.Arrays;
import java.util.ArrayList;
/**
 * Write a description of class IFFT here.
 *
 * @author (your name)
 * @version (a version number or a date)
 */
public class IFFT
{

    public IFFT()
    {

    }

    public static void test(){
        float[] sound = ReadSound.readSound("ss.wav");
        double[] sig = new double[sound.length];
        for(int i = 0; i < sig.length; i++){
            sig[i] = sound[i];
        }
        ArrayList<double[]> spect = FFT.getSpectrumWPhase(Arrays.copyOf(sig, WaveWriter.SAMPLE_RATE * 3), false);
        double[] fDomain = new double[spect.size()];
        for(int i = 0; i < fDomain.length; i++){
            fDomain[i] = spect.get(i)[1];
        }
        spect = FFT.getSpectrumWPhase(fDomain, false);
        sig = new double[fDomain.length];
        for(int i = 0; i < sig.length; i++){
            sig[i] = spect.get(i)[1];
        }
        WaveWriter ww = new WaveWriter("ifft");
        sig = normalize(sig);
        for(int i = 0; i < sig.length; i++){
            ww.df[0][i] += sig[i];
            ww.df[1][i] += sig[i];
        }
        ww.render();
    }

    private static double[] normalize(double[] sig){
        double max = Double.MIN_VALUE;
        double min = Double.MAX_VALUE;
        for(int i = 0; i < sig.length; i++){
            max = Math.max(max, sig[i]);
            min = Math.min(min, sig[i]);
        }

        for(int i = 0; i < sig.length; i++){
            sig[i] = 2 * (sig[i] - min) / (max - min)  - 1;
        }
        return sig;
    }
}
