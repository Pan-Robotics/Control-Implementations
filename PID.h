#ifndef PID_H
#define PID_H

class PID
{
private:
    float gains[3];
    float tau;
    float dt;
    float error;
    float sumError;
    float prevError;
    float prevY;
    float prevU;

public:
    // Constructor
    PID(float Kp, float Ki, float Kd, float Tau)
    {
        gains[0] = Kp;
        gains[1] = Ki;
        gains[2] = Kd;
        tau = Tau;
        dt = 0.0f;
        error = 0.0f;
        sumError = 0.0f;
        prevError = 0.0f;
        prevY = 0.0f;
        prevU = 0.0f;
    }

    // Simple PID With Derivative Filter
    float pid(float Val, float Sp)
    {
        error = Sp - Val;

        float P = gains[0] * error;

        float I = gains[1] * sumError;
        sumError += dt * error;

        float D = gains[2] * filter(error - prevError) / dt;

        prevError = error;

        float sumPID = saturation(P + I + D, -1000, 1000);
        return sumPID;
    }

    // Simple Bilinear First Order Filter
    float filter(float U)
    {
        float coeff = 2 * tau / dt;
        float num = 1 - coeff;
        float dnum = 1 + coeff;
        float Y = (U + prevU - prevY * num) / dnum;

        prevU = U;
        prevY = Y;
        return Y;
    }

    // 'delta' is Elapsed Time Since Previous Frame.
    void pid_run(float delta)
    {
        dt = delta;
    }

    float saturation(float data, float min_val, float max_val)
    {
        if (data > max_val)
        {
            data = max_val;
        }
        if (data < min_val)
        {
            data = min_val;
        }
        return data;
    }
};

#endif // PID_H
