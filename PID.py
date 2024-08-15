class PID:
    def __init__(self, Kp, Ki, Kd, Tau):
        self.gains = [Kp, Ki, Kd]
        self.tau = Tau
        self.dt = 0.0
        self.error = 0.0
        self.sum_error = 0.0
        self.prev_error = 0.0
        self.prev_y = 0.0
        self.prev_u = 0.0

    def pid(self, val, sp):
        self.error = sp - val

        P = self.gains[0] * self.error

        I = self.gains[1] * self.sum_error
        self.sum_error += self.dt * self.error

        D = self.gains[2] * self._filter(self.error - self.prev_error) / self.dt

        self.prev_error = self.error

        sum_pid = self._saturation(P + I + D, -1000, 1000)
        return sum_pid

    def _saturation(self, data, min_val, max_val):
        if data > max_val:
            data = max_val
        if data < min_val:
            data = min_val
        return data

    def _filter(self, u):
        coeff = 2 * self.tau / self.dt
        num = 1 - coeff
        dnum = 1 + coeff
        y = (u + self.prev_u - self.prev_y * num) / dnum

        self.prev_u = u
        self.prev_y = y
        return y

    def pid_run(self, delta):
        self.dt = delta
