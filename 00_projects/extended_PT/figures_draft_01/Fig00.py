from back_process import *


if __name__ == "__main__":

    x = np.arange(0, 2 * np.pi, 0.01)
    #x = np.arange(0, 360, 0.01)
    y = np.empty_like(x)

    m_1 = 0 #25 / np.pi
    m_2 = - 12 / np.pi
    m_3 = 12 / np.pi

    gtr = np.pi / 180
    zero = 28

    for i in range(len(x)):
        if 0 * gtr <= x[i] < 45 * gtr:
            y[i] = m_1 * x[i] + zero
        elif 45 * gtr <= x[i] < 80 * gtr:
            y[i] = float("NaN")
        elif 80 * gtr <= x[i] < 165 * gtr:
            y[i] = m_2 * x[i] - m_2 * 120 * gtr + zero
        elif 165 * gtr <= x[i] < 195 * gtr:
            y[i] = float("NaN")
        elif 195 * gtr <= x[i] < 285 * gtr:
            y[i] = m_3 * x[i] - m_3 * 240 * gtr + zero
        elif 285 * gtr <= x[i] < 310 * gtr:
            y[i] = float("NaN")
        elif 310 * gtr <= x[i] < 360 * gtr:
            y[i] = m_1 * x[i] - m_1 * 360 * gtr + zero
    df = pd.DataFrame(y, columns=['a'])
    df = df.interpolate(method='polynomial', order=5)

    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.plot(x, df, color="r", lw=2.5)
    ax.fill_between(np.arange(75 * gtr, 165 * gtr, 0.25), 0, 50, color="r", alpha = 0.2)
    ax.fill_between(np.arange(195 * gtr, 285 * gtr, 0.25), 0, 50, color="r", alpha=0.2)
    ax.set_rmax(35)
    ax.set_rticks([10, 20, 30])# Less radial ticks
    ax.tick_params(labelsize=15)
    ax.set_rlabel_position(0)  # Move radial labels away from plotted line
    plt.show()