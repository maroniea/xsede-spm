from dataclasses import dataclass, asdict
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import base64
import re
from scipy import optimize
from munch import Munch
from util import process_file

import lmfit as lm
import sympy as sm
import schemdraw
import schemdraw.elements as elm

import io


def base64png(fig):
    pic_IObytes = io.BytesIO()
    fig.savefig(pic_IObytes, format='png')
    pic_IObytes.seek(0)
    pic_hash = base64.b64encode(pic_IObytes.read())
    return pic_hash.decode('utf-8')


def base64svgfig(fig, alt):
    pic_IObytes = io.BytesIO()
    fig.savefig(pic_IObytes, format='svg')
    pic_IObytes.seek(0)
    b64 = base64.b64encode(pic_IObytes.read()).decode('utf-8')
    return f'<img src="data:image/svg+xml;base64, {b64}" alt="{alt}"/>'


def base64svg(svg, alt):
    b64 = base64.b64encode(svg).decode("utf-8")
    return f'<img src="data:image/svg+xml;base64, {b64}" alt="{alt}"/>'

def img64(b64, alt="Report figure"):
    return f'<img src="data:image/png;base64, {b64}" alt="{alt}"/>'


def html_output(fig, f2, f_r, d2, fits, files, data, cis):
    fig_png = base64png(fig)
    f2_png = base64png(f2)
    fname = files[0].name
    report = fits[0].fit_report()
    ci_report = lm.ci_report(cis[0])
    raw_data = data[0].to_html(classes='mystyle')
    html = f"""
<html>
<head>
<style>
.mystyle {{
font-size: 11pt; 
font-family: Arial;
border-collapse: collapse; 
border: 1px solid silver;

}}

.mystyle td, th {{
padding: 5px;
}}

.mystyle tr:nth-child(even) {{
background: #E0E0E0;
}}

.mystyle tr:hover {{
background: silver;
cursor: pointer;
}}
</style>
</head>
<body>
<h1>Streamlit Impedance Report</h1>

<p>Filename: {fname}</p>

<h3>Parameters and Fit Report</h3>
<pre>
<code>
{report}
</code>
</pre>

<h4>Confidence Intervals</h4>
<pre>
<code>
{ci_report}
</code>
</pre>

<h2>Circuit</h2>
{base64svg(d2.get_imagedata(), "Circuit")}


<h2>Impedance vs. Frequency</h2>
{base64svgfig(fig, alt="Impedance")}

<h3>Residuals</h3>
{base64svgfig(f_r, alt="Impedance Resid.")}

<h2>Nyquist Plot</h2>
{base64svgfig(f2, alt="Impedance")}




<h3>Raw Data</h3>

{raw_data}
</body>
</html>"""
    return html

    


# import plotly.express as px
# import pandas as pd


class ZElements:
    def R(f, R):
        return R
    
    def C(f, C):
        omega = 2*np.pi*f
        return 1/(1j*omega*C)
    
    def W(f, Aw):
        sqrt_omega = np.sqrt(2*np.pi*f)
        return Aw/sqrt_omega + Aw/(1j*sqrt_omega)




Zdef = dict(R=100.0, C=1.0, W=10.0, P=1.0, n=0.5, L=10.0)

def impedance(x):
    try:
        return Z[x.name[0]](x)
    except:
        return x

def default_val(x):
    return Zdef[x[0]]



class Model:
    
    def __init__(self):
        self.params = set()

    def par(self, x, y):
        return 1/(1/self.impedance(x) + 1/self.impedance(y))

    def ser(self, *args):
        return sum(self.impedance(x) for x in args)

    def impedance(self, x):
        try:
            imp = Z[x.name[0]](x)
            if x not in self.params:
                self.params.append(x)
        except:
            imp = x
        
        return imp

def render_svg(svg):
    """Renders the given svg string."""
    b64 = base64.b64encode(svg).decode("utf-8")
    html = r'<img src="data:image/svg+xml;base64,%s"/>' % b64
    st.write(html, unsafe_allow_html=True)

cDict = dict(R=elm.Resistor, C=elm.Capacitor, W=elm.RBox, P=elm.RBox, L=elm.Inductor)



class EClass:
    f = sm.symbols('f')
    Z = dict(R=lambda x: x, C=lambda x: 1/(2*sm.pi*sm.I*EClass.f*x),
        W= lambda x: x/sm.sqrt(2*sm.pi*EClass.f)+x/(1j*sm.sqrt(2*sm.pi*EClass.f)),
        P= lambda x, y: 1/(x*(2*sm.pi*EClass.f)**y)*sm.exp(-sm.pi/2*y*sm.I),
        L=lambda x: 2*sm.I*sm.pi*EClass.f*x
        )

    def __call__(self, label):
        if label[0] != "P":
            symbol = sm.symbols(label)
            Z_val = self.Z[label[0]](symbol)
            im = ImpedanceMod(Z_val, params={symbol}, network=label)
        else:
            P, n = sm.symbols(f"{label} n_{label[1:]}")
            Z_val = self.Z[label[0]](P, n)
            im = ImpedanceMod(Z_val, params={P, n}, network=label)
        
        return im


E = EClass()

# When we have sub, we can just append...

# When we have par, we need to add a wire up and down, then put



@dataclass
class Parallel:
    data: list

    def draw(self, d, space=0.4):
        first, second = self.data
        d.push()
        xi, yi = d.here

        delta_len = abs(len(second) - len(first))
        # print(delta_len)
        if self.first_branch_longer():
            d += elm.Line().up(d.unit*space)
            d += elm.Line().right(d.unit*0.01)
            draw_string(first, d)
            d += (final_pos := elm.Line().down(d.unit*space))
            xf1, yf1 = d.here
            # Back to previous position
            d.pop()
            tot = (xf1 - xi)
            d += elm.Line().down(d.unit*space)
            d += elm.Line().right(d.unit*0.01+delta_len*0.5*tot)
            draw_string(second, d)
        else:
            d += elm.Line().down(d.unit*space)
            d += elm.Line().right(d.unit*0.01)
            draw_string(second, d)
            d += (final_pos := elm.Line().up(d.unit*space))
            xf1, yf1 = d.here
            d.pop()
            tot = (xf1 - xi)
            d += elm.Line().up(d.unit*space)
    
            d += elm.Line().right(d.unit*0.01+d.unit*delta_len*0.5)
            draw_string(first, d)

        
        d += elm.Wire('-|').to(final_pos.end)           
        d += elm.Line().right(d.unit*space)


    def __len__(self):
        left, right = self.data
        left_len = 1 if isinstance(left, str) else len(left)
        right_len = 1 if isinstance(right, str) else len(right)
        return max(left_len, right_len)
    
    def first_branch_longer(self):
        left, right = self.data
        left_len = 1 if isinstance(left, str) else len(left)
        right_len = 1 if isinstance(right, str) else len(right)
        return left_len > right_len


def parallel(*args):
    return Parallel(args)

@dataclass
class Series:
    data: list

    def draw(self, d, space=0.4):
        
        for i, x in enumerate(self.data):
            draw_string(x, d)
    

    def __len__(self):
        s = 0
        for x in self.data:
            if isinstance(x, str):
                s+=1
            else:
                s+=len(x)
        return s


def series(*args):
    return Series(args)

def draw_string(x, d):
    if isinstance(x, (Parallel, Series)):
        x.draw(d)
    else:
        element = cDict[x[0]]
        d += element(label=x)

class ImpedanceMod:

    f = sm.symbols('f')

    def __init__(self, Z, params, network):
        self.Z = Z
        self.params = set(params)
        self.network = network

        self.Zfunc = sm.lambdify([self.f, *self.params_sorted],
                        self.Z.subs(self.f, self.f*1e-6), 'numpy')



    def create_model(self, **kwargs):
        self.ZModel = lm.Model(self.Zfit)
        for p in self.param_names:
            if p in kwargs:
                val = kwargs[p]
            elif p[0] == 'n':
                val = 0.5
            else:
                val = 1.0
            if p[0] == 'n':
                self.ZModel.set_param_hint(p, value=val, vary=True, min=0.0, max=1.0)
            else:
                self.ZModel.set_param_hint(p, value=val, vary=True, min=0.0)

        self.pars = self.ZModel.make_params()
        return self.pars

    def Zfit(self, f, **params):
        # Return real, then imaginary parts...
        Z = self.Zkw(f, **params)
        return np.r_[Z.real, Z.imag]
    
    def Zkw(self, f, **params):
        return self.Zfunc(f, *[params[key] for key in self.param_names])
        

    def __floordiv__(self, other):
        Z_new = 1/(1/self.Z + 1/other.Z)
        return ImpedanceMod(Z_new, self.params|other.params,
                            parallel(self.network,other.network)
                            )
    
    def __sub__(self, other):
        Z_new = self.Z + other.Z
        return ImpedanceMod(Z_new, self.params|other.params,
                            series(self.network, other.network))
    
    def __repr__(self):
        return 'ImpedanceMod('+repr(self.Z)+','+repr(self.params)+')'
    
    @property
    def params_sorted(self):
        return sorted(list(self.params), key= lambda x: x.name)

    @property
    def param_names(self):
        return [p.name for p in self.params_sorted]
    
def circuit_from_string(string):
    str_to_eval = re.sub('([a-zA-Z0-9]+)', "E('\\1')", string)
    return eval(str_to_eval)


def run():

    if 'ever_submitted' not in st.session_state:
        st.session_state.ever_submitted = False

    st.title("Fit Impedance Model")

    st.markdown("""
This page fits electrochemical impedance spectroscopy data to a circuit chosen below.
""")

    st.markdown("## Circuit")

    # Here is a nice circuit - we should probably be saving the value to local storage

    circuit_str = st.text_input("Write out the circuit:", value="R1-R2//C2")

    with st.expander("Circuit help?"):
        st.markdown("""To use a circuit element, choose the first letter of the
        string to be
        
        R: Resistor (Ω)
        C: Capacitor (uF)
        L: Inductor (uH)
        P: Constant Phase Element
            - P (commonly Q, units Ω⁻¹ μsⁿ)
            - n: Unitless, 0 < n < 1
        W: Warburg Element (Ω)
        
To connect circuit elements, use `-` for a series connection and
`//` for a parallel connection. You can use parenthesis to clarify
any ambiguities. For example,
R1 and R2 (in series) in parallel with
C1 would be written `(R1-R2)//C1`""")

    mod = circuit_from_string(circuit_str)


    mod.create_model()

    with schemdraw.Drawing(show=False) as d2:
        mod.network.draw(d2)

    render_svg(d2.get_imagedata())
    st.markdown("Impedance:")
    st.write(mod.Z)
        

    st.markdown("""## Load Experimental Data""")

    files = [st.file_uploader("Load Experimental Data", accept_multiple_files=False)]

    data = []
    labels = []
    # Process each file...

    if None not in files:
        st.write(files)

        filenames = [(i, f.name) for i, f in enumerate(files)]
        data = [process_file(f) for f in files]

        ind_fname = st.selectbox("Choose data to display: ", filenames,
            format_func=lambda x: x[1], index=0)

        st.write("""### Labels
Use the boxes below to change the labels for each line that will go on the graph.""")
        labels = [st.text_input(f"{filename[0]}. {filename[1]}", value="") for filename in filenames]
        
        if ind_fname:
            df = data[ind_fname[0]]
            cols = list(df.columns)


        st.write("### Choose columns")
        with st.form("column_chooser_and_run"):
            f_column = st.selectbox("Choose the frequency column: ", cols)
            x_column = st.selectbox("Choose Z' (in-phase, X) column: ", cols, index=len(cols)-2)
            y_column = st.selectbox("Choose Z'' (out-of-phase, Y) column: ", cols, index=len(cols)-1)

            submitted = st.form_submit_button()

        
        st.session_state.ever_submitted = submitted | st.session_state.ever_submitted

    
    plotModel = st.sidebar.checkbox("Plot initial model?", value=True)
    if plotModel:
        st.sidebar.markdown("### Frequencies")
        freqs = np.geomspace(
                st.sidebar.number_input("Start freq (Hz)", value=1.0, format="%.1e"),
                st.sidebar.number_input("End Freq (Hz)", value=1e4, format="%.1e"),
                int(st.sidebar.number_input("Pts", value=100, format="%d"))
        )

    st.sidebar.markdown("### Initial Parameters")
    param_values = [st.sidebar.number_input(x.name, value=default_val(x.name),
                        format='%.2e') for x in mod.params_sorted]

    if plotModel:
        out = mod.Zfunc(freqs, *param_values)

    fits = []
    cis = []
    for label, d in zip(labels, data):
        y_data = np.r_[d[x_column].values, d[y_column].values]
        param_items = dict(zip(mod.param_names, param_values))
        pars = mod.create_model(**param_items)
        f_data = d[f_column].values
        result = mod.ZModel.fit(y_data, pars, f=f_data)
        def func(kwargs):
            return mod.ZModel.eval(f=f_data, **kwargs) - y_data
        minimizer = lm.Minimizer(func, result.params)
        ci = lm.conf_interval(minimizer, result, sigmas=[1, 2])
        cis.append(ci)
        fits.append(
            result
        )
    
    st.markdown("""## Results""")

    for fit, ci in zip(fits,cis):
        st.text(fit.fit_report())
        ci_report = lm.ci_report(ci)
        st.text(ci_report)

    f2, a2 = plt.subplots()

    for d, fit, label in zip(data, fits, labels):
        # st.write(fit.best_values)
        Z_fit = mod.Zfunc(d[f_column].values, **fit.best_values)
        a2.plot(d[x_column], -d[y_column], '.', label=label)
        a2.plot(Z_fit.real, -Z_fit.imag , label=label+' fit')
    
    if plotModel:
        a2.plot(out.real, -out.imag, label='Model', color='0')

    a2.set_xlabel("$ Z' $ (ohm)")
    a2.set_ylabel("$-Z''$ (ohm)")
    a2.set_aspect('equal', 'box')
    a2.legend()
    f2.tight_layout()

    st.pyplot(f2)


    fig, ax = plt.subplots()

    f_r, ax_r = plt.subplots()

    color_cycle = ax._get_lines.prop_cycler
    for d, fit, label in zip(data, fits, labels):
        Z_fit = mod.Zfunc(d[f_column].values, **fit.best_values)
        l1, = ax.semilogx(d[f_column], d[x_column], '.', label=label+' Z\'')
        new_color = next(color_cycle)['color']
        l2 = ax.scatter(d[f_column], -d[y_column], s=6, marker='o',fc='none', ec=new_color, label=label+' -Z\"')
        
        ax.semilogx(d[f_column], Z_fit.real, '-', color=l1.get_color(), label=label+" $Z'$ mod")
        ax.semilogx(d[f_column], -Z_fit.imag, '--', color=new_color, label=label+" $-Z''$ mod")

        ax_r.semilogx(d[f_column], d[x_column].values - Z_fit.real, '.', color=l1.get_color(), label="$Z'$")
        ax_r.scatter(d[f_column], -d[y_column].values + Z_fit.imag, s=6, marker='o', ec=new_color, fc='none', label="$Z''$")

    # Model
    if plotModel:
        ax.semilogx(freqs, out.real, color='0', label="Model $Z'$")
        ax.semilogx(freqs, -out.imag,  '--', color='0',label="Model $-Z''$")
    
    for ax_ in [ax, ax_r]:
        ax_.set_xlabel("Frequency (Hz)")
        ax_.legend()

    ax.set_ylabel("$Z$ (ohm)")
    ax_r.set_ylabel("Resid. $r = Z - Z_\\mathrm{model}$(ohm)")

    st.pyplot(fig)

    st.pyplot(f_r)

    if None not in files:
        st.download_button(label="Download Report", data=html_output(fig, f2, f_r, d2, fits, files, data, cis), file_name='report.html')
    

if __name__ == '__main__':
    run()