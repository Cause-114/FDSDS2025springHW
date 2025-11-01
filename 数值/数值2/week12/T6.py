import numpy as np
import matplotlib.pyplot as plt
import sounddevice as sd
import soundfile as sf
import threading
import time
import queue

# 这个地方我们只能取一个通道（一般的音频文件都是多通道）
# 但是老师那个电话铃声示例就是单通道，不用处理
# y, sr = sf.read("./week11/DTMF_dialing.ogg")
y, sr = sf.read("114.mp3")
y = y[:, 0]  # 取第一列（左通道）

# 每一次频谱分析需要的数据点个数
frame_size = 2048

# 创建图形窗口
plt.ion()  # 启用交互模式
fig, ax = plt.subplots()
line, = ax.plot([], [], 'b')
text = ax.text(0.5, 0.9, '', transform=ax.transAxes, ha='center')
# 这里的line跟text相当于占位符，后面会被频谱数据替换
ax.set_ylim(0, 100)
ax.set_xlim(0, sr // 2)
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Magnitude')

# 创建线程安全的队列
data_queue = queue.Queue()

# 频谱分析函数
def compute_fundamental(signal, sr):
    # 这个window据说是对首尾“磨光”，保证效果更准。
    window = np.hanning(len(signal))
    # np.fft.rfft只会输出正频率（由于实数的DFT转换到到
    # 三角插值时正负频率系数共轭，所以只用输出一半）
    spectrum = np.abs(np.fft.rfft(signal * window))
    # 输出与上面的一半DFT相匹配的频率
    freqs = np.fft.rfftfreq(len(signal), d=1/sr)
    # 找到最大值对应的索引
    peak_index = np.argmax(spectrum)
    fundamental = freqs[peak_index]
    return fundamental, freqs, spectrum

pos = 0
play_done = False  # 添加标志位，控制音频播放结束

def callback(outdata, frames, time, status):
    global pos, play_done
    chunk = y[pos:pos+frames]
    # chunk代表当前窗口的音频数据（要输出到播放端的）
    if len(chunk) < frames:
        outdata[:len(chunk), 0] = chunk
        outdata[len(chunk):, 0] = 0
        play_done = True  # 音频播放完毕
        raise sd.CallbackStop()  # 停止回调
    else:
        outdata[:, 0] = chunk

    # 截取输出数据中的一部分做频谱分析
    segment = chunk[:frame_size]
    if len(segment) > 0:  # 只要segment有数据，就进行频谱分析
        f0, freqs, spectrum = compute_fundamental(segment, sr)
        # 将数据放入队列
        data_queue.put((freqs, spectrum, f0))

    pos += frames

def play_audio():
    # 启动播放，由系统调用callback函数，outdata是系统想要知道的输出信息，frames是系统想要的样本点个数
    with sd.OutputStream(callback=callback, samplerate=sr, channels=1, dtype='float32'):
        while not play_done:  # 检查播放是否完成
            time.sleep(0.01)  # 等待 0.1 秒
# 启动音频播放线程
audio_thread = threading.Thread(target=play_audio)
audio_thread.start()

# 图形更新
while not play_done:
    if not data_queue.empty():
        freqs, spectrum, f0 = data_queue.get()  # 从队列中获取频谱数据
        line.set_data(freqs, spectrum)  # 更新频谱图
        text.set_text(f"Fundamental: {f0:.2f} Hz")  # 显示基本频率
        fig.canvas.draw()  # 绘制
        fig.canvas.flush_events()  # 刷新事件循环
    
    plt.pause(0.01)  # 每 0.01 秒更新一次图形

# 确保图形在音频播放结束后保持显示
plt.ioff()  # 关闭交互模式，防止图形在播放时消失
plt.show()  # 显示图形
