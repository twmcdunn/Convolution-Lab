import java.io.InputStream;
import java.io.IOException;
public class StereoPcmInputStream extends InputStream
{
    private float[][] dataFrames;
    private int framesCounter;
    private int cursor;
    private int[] pcmOut = new int[2];
    private int[] frameBytes = new int[4];
    private int idx;

    public int framesToRead;

    public void setDataFrames(float[][] dataFrames)
    {
        this.dataFrames = dataFrames;
        framesToRead = dataFrames[0].length;// / 2;
    }

    @Override
    public int read() throws IOException
    {
        while(available() > 0)
        {
            idx &= 3; 
            if (idx == 0) // set up next frame's worth of data
            {
                framesCounter++; // count elapsing frames

                // scale to 16 bits
                pcmOut[0] = (int)(dataFrames[0][cursor] * Short.MAX_VALUE);
                pcmOut[1] = (int)(dataFrames[1][cursor++] * Short.MAX_VALUE);

                // output as unsigned bytes, in range [0..255]
                frameBytes[0] = (char)pcmOut[0];
                frameBytes[1] = (char)(pcmOut[0] >> 8);
                frameBytes[2] = (char)pcmOut[1];
                frameBytes[3] = (char)(pcmOut[1] >> 8);

            }
            return frameBytes[idx++]; 
        }
        return -1;
    }

    @Override 
    public int available()
    {
        // NOTE: not concurrency safe.
        // 1st half of sum: there are 4 reads available per frame to be read
        // 2nd half of sum: the # of bytes of the current frame that remain to be read
        return 4 * ((framesToRead - 1) - framesCounter) 
                + (4 - (idx % 4));
    }    

    @Override
    public void reset()
    {
        cursor = 0;
        framesCounter = 0;
        idx = 0;
    }

    @Override
    public void close()
    {
        System.out.println(
            "StereoPcmInputStream stopped after reading frames:" 
                + framesCounter);
    }
}