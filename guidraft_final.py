from tkinter import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
root=Tk()
root.title("H11 graph calculator")
# making the labels 
label1=Label(root,text="radius of tread ring (m)")
label1.grid(row=0,column=0,columnspan=1)

label2=Label(root,text="tread width (m)")
label2.grid(row=1,column=0,columnspan=1)

label3=Label(root,text="density of tread ring (kg/m)")
label3.grid(row=2,column=0,columnspan=1)

label4=Label(root,text="mass of wheel (kg)")
label4.grid(row=3,column=0,columnspan=1)

label5=Label(root,text=" moment of inertia of wheel (kg·m²)")
label5.grid(row=4,column=0,columnspan=1)

label6=Label(root,text="flexural rigidity of tread ring (N·m²)")
label6.grid(row=5,column=0,columnspan=1)

label7=Label(root,text="radial fundamental spring rate (N/m²)")
label7.grid(row=6,column=0,columnspan=1)

label8=Label(root,text="circumferential fundamental spring rate (N/m²)")
label8.grid(row=7,column=0,columnspan=1)

label9=Label(root,text="inflation pressure (Pa)")
label9.grid(row=8,column=0,columnspan=1)

label10=Label(root,text="circumferential tension of tread ring (N)")
label10.grid(row=9,column=0,columnspan=1)

label11=Label(root,text="damping coefficient")
label11.grid(row=10,column=0,columnspan=1)

label12=Label(root,text="angular velocity (rad/s)")
label12.grid(row=11,column=0,columnspan=1)

# making entry boxes for the values 
e1=Entry(root)
e1.grid(row=0,column=1,columnspan=3)

e2=Entry(root)
e2.grid(row=1,column=1,columnspan=3)

e3=Entry(root)
e3.grid(row=2,column=1,columnspan=3)

e4=Entry(root)
e4.grid(row=3,column=1,columnspan=3)

e5=Entry(root)
e5.grid(row=4,column=1,columnspan=3)

e6=Entry(root)
e6.grid(row=5,column=1,columnspan=3)

e7=Entry(root)
e7.grid(row=6,column=1,columnspan=3)

e8=Entry(root)
e8.grid(row=7,column=1,columnspan=3)

e9=Entry(root)
e9.grid(row=8,column=1,columnspan=3)

e10=Entry(root)
e10.grid(row=9,column=1,columnspan=3)

e11=Entry(root)
e11.grid(row=10,column=1,columnspan=3)

e12=Entry(root)
e12.grid(row=11,column=1,columnspan=3)
def graph_maker():
    try:
        a = float(e1.get())
        b = float(e2.get())
        rho_A = float(e3.get())
        m_wheel = float(e4.get())
        I_r = float(e5.get())
        EI = float(e6.get())
        k_r = float(e7.get())
        k_t = float(e8.get())
        p_0 = float(e9.get())
        sigma_theta_0_A = float(e10.get())
        zeta = float(e11.get())
        Omega = float(e12.get())

        def t_m0_11(s):
            m_0 = rho_A
            k_0R = k_t * a
            k_R = k_t * a**2
            m_r = I_r / (2 * np.pi * a)
            c_0 = 2 * zeta * np.sqrt(m_0 * k_0R)
            c_r = 2 * zeta * np.sqrt(m_r * k_R)
            num = m_r * s**2 + c_r * s + k_R
            den = (m_0 * s**2 + c_0 * s + k_0R) * (m_r * s**2 + c_r * s + k_R) - k_0R * 2
            return num / den

        def t_m1_11(s):
            m_1 = rho_A * 2
            m_a = m_wheel / (np.pi * a)
            k_1 = k_r + k_t - rho_A * Omega**2 * 2
            k_a = k_r + k_t - (m_wheel * Omega**2) / (np.pi * a)
            c_1 = 2 * zeta * np.sqrt(m_1 * k_1)
            c_a = 2 * zeta * np.sqrt(m_a * k_1)
            num = m_a * s**2 + c_a * s + k_a
            den = (m_1 * s**2 + c_1 * s + k_1) * (m_a * s**2 + c_a * s + k_a) - k_1**2
            return num / den

        def t_mn_11(n, s):
            if n == 0:
                return t_m0_11(s)
            elif n == 1:
                return t_m1_11(s)
            else:
                m_n = rho_A * (1 + n**2)
                k_n = ((EI / a**4) * n**2 + (sigma_theta_0_A / a**2)) * (n**2 - 1)**2 + (p_0 * b / a) * (n**2 - 1) + k_r * n**2 + k_t - rho_A * Omega**2 * (n**2 + 1)
                c_n = 2 * zeta * np.sqrt(m_n * k_n)
                g_n = 2 * rho_A * n * Omega * (n**2 - 1)
                num = m_n * s**2 + c_n * s + k_n
                den = (m_n * s**2 + c_n * s + k_n)**2 + (g_n * s)**2
                return num / den

        def h11(omega, n_max=10):
            s = 1j * omega
            return (0.5 * t_m0_11(s) + sum(t_mn_11(n, s) for n in range(1, n_max + 1))) / (np.pi * a)

        
        frequencies = np.linspace(0.1, 300, 200)
        h11_magnitudes = [abs(h11(2 * np.pi * f)) for f in frequencies]

        
        global graph_canvas
        if 'graph_canvas' in globals():
            graph_canvas.get_tk_widget().destroy()

        # Create figure
        fig = Figure(figsize=(6, 4), dpi=100)
        ax = fig.add_subplot(111)
        ax.loglog(frequencies, h11_magnitudes)
        ax.set_title('H11 Transfer Function Magnitude')
        ax.set_xlabel('Frequency (Hz)')
        ax.set_ylabel('|H11| Magnitude')
        ax.grid(True, which="both", ls="--")

        # Embed it
        graph_canvas = FigureCanvasTkAgg(fig, master=root)
        graph_canvas.draw()
        graph_canvas.get_tk_widget().grid(row=0, column=4, rowspan=15, padx=10, pady=10)

    except ValueError:
        from tkinter import messagebox
        messagebox.showerror("Input Error", "Please enter valid numeric values for all parameters.")

def clear():
    e1.delete(0,END)
    e2.delete(0,END)
    e3.delete(0,END)
    e4.delete(0,END)
    e5.delete(0,END)
    e6.delete(0,END)
    e7.delete(0,END)
    e8.delete(0,END)
    e9.delete(0,END)
    e10.delete(0,END)
    e11.delete(0,END)
    e12.delete(0,END)
    return
button1=Button(root,text="Graph It!!",command=graph_maker)
button1.grid(row=12,column=5,columnspan=2)

button2=Button(root,text="Clear All",command= clear)
button2.grid(row=13,column=5,columnspan=2)
root.mainloop()