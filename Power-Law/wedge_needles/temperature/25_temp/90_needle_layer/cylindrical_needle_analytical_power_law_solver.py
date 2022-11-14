import numpy              as np
import pandas             as pd
import matplotlib.pyplot  as plt
import matplotlib.patches as patches
import matplotlib.ticker  as ticker
from   matplotlib         import cm
from   math               import pi
from   scipy              import integrate, stats
from   scipy.misc         import derivative

# mathematical plot font
#import matplotlib
#matplotlib.rcParams['mathtext.fontset'] = 'stix'
#matplotlib.rcParams['font.family'] = 'STIXGeneral'

version_number = "3.2.6"
update_content ="corrected printing output"
print('################################################################################################')
print( f' Bioink inside a cylindrical needle; analytical solver based on power law, KORIN, version {version_number}\n')
print( ' #UPDATE LOG#\n', f'-{update_content}')
print('################################################################################################\n')

# DIY SECTION #
###########################
# needle parameters

# needle radius (input: micro-m) [unit: m]
R   = 400 * (10**-6)

# needle length (input: mm) [unit: m]
Ln  = 13 * (10**-3)

# pressure or flow rate is known
pressure_is_known = False
if pressure_is_known == True:
    # pressure drop along the needle (input: kPa) [unit: Pa]
    Pn  = 200 * (10**3)
else:
    # volumetric flow rate (input: micro-L/s) [unit: m^3/s]
    Qn  = 50 * (10**-9)

###########################
# bioink properties (power law)

# constant flow index [Pa*s^n]
K   = 29.86

# flow behavior index [-]
n   = 0.46

# percent weight/volume (% w/v)
#wvp = 2

###########################

# paraview csv file name (place file inside the same folder)
# Paraview -> slice over needle cross section -> plot over line -> set point 0 to be on the center -> save data as .csv
data_file_name = "90_data_25C" # do not put filename extension
plot_graphs = True
save_graphs = False
dpi_save = 100;
#############################################################################################

# CONSTANT SECTION #

#rho_solution = 998.23              # solution desnsity, in this case water@20C [kg/m^3]
#rho = rho_solution + wvp*10        # bioink density if the solution is water@20C [kg/m^3]
rho = 1000

#############################################################################################

# EQUATION SECTION #

# velocity profile along the needle variable radius, Vz(r):
def Vz(r):
    return ((n/(n+1))*((Pn*R)/(2*K*Ln))**(1/n)*R*(1-(r/R)**((n+1)/n)))

# product of r and Vz(r), r*Vz(r):
#def rVz(r):
#	return r*Vz(r)

# average volumetric flow rate, Q(R):
def Qave(R):
    return (8*pi*R**3)/(8)*((Pn*R)/(2*K*Ln))**(1/n)*(n/(3*n+1))

# pressure drop based on the average volumetric flow rate, Pn
def Pn_func(Qn):
    return (Qn/((8*pi*R**3)/(8)*((R)/(2*K*Ln))**(1/n)*(n/(3*n+1))))**n

#############################################################################################

# CALCULATION SECTION #

N = 1000 # linspace step
x = np.linspace(1e-6, R, N) # numpy array from 0 to needle radius, R

"""
# find flow rate, Q(r),: # without initial condition Vz_0
Q = []
x_range = []
y_range = []
for i in x:
    x_range.append(i)
    y_range.append(rVz(i))
    integral = integrate.simps(y_range, x_range) # integral of 2*pi*r*Vz(r) dr from r to R
    Q.append(integral)
Q = list(map(lambda var: var*2*pi, Q)) # times conatants (2 and pi) to get Q(r)
"""
if pressure_is_known == True:
    Pn = Pn
else:
    Pn = Pn_func(Qn) # find pressure drop from flow rate

tol = 1e-6 # tolerance
dVzdr = derivative(Vz, x, tol) # derivative of Vz(r) with respect to r

# POWER LAW #

def power_n(my_negative_list, n_minus_1 = 0): # raise power of each items inside a list
    if n_minus_1 == -1:
        list_postive = [j**(n-1) for j in -1*my_negative_list] # avoid rasing power for negative numbers
        return list(map(lambda var: var*-1, list_postive)) # change the signs back again
    else:
        list_postive = [j**n for j in -1*my_negative_list] # avoid rasing power for negative numbers
        return list(map(lambda var: var*-1, list_postive)) # change the signs back again

gamma_dot_n = power_n(dVzdr) # shear rate raised to the power of n
tau_rz = list(map(lambda var: var*-K, gamma_dot_n)) # obtain shear stress (absoulte value)

way_power_law_option = 2 # input = 1 or 2 (same result)

if way_power_law_option == 1:
    # way 1 to obtain viscosity
    gamma_dot_n_minus_1 = power_n(dVzdr, -1) # shear rate raised to the power of n-1
    eta = list(map(lambda var: var*-K, gamma_dot_n_minus_1)) # obtain apparent viscosity
if way_power_law_option == 2:
    # way 2 to obtain viscosity
    eta = tau_rz/-dVzdr # obtain apparent viscosity (make dVzdr positive since tau_rz is now positive)

#############################################################################################

# PRINTING SECTION #

Q_average    = Qave(R)
Q_average_mL = Qave(R)*1e6
print("Volumetric Flow Rate =", Q_average_mL, "[mL/s]\n")

Q_ave_mass_mg = Q_average*rho*1e3
print("Mass Flow Rate =" , Q_ave_mass_mg, "[mg/s]\n")

Pn_kPa = Pn/1e3
print("Pressure Drop along the Needle =", Pn_kPa, "[kPa]\n")

Vz_max_mm = Vz(0)*1e3
V_average = (Pn/2/K/Ln)**(1/n)*(n/(3*n+1))*R**((n+1)/n) # Qn/pi/R**2
print("Needle Center Line Velocity =",Vz_max_mm, "[mm/s]; Average Extrusion Velocity =",V_average*1e3,"[mm/s]\n")

tau_max_kPa = (R/2)*(Pn/Ln)/1e3
print("Wall Shear Stress =", tau_max_kPa, "[kPa]\n")

print("Max Shear Rate =", max(-dVzdr), "[1/s];" , "Min Shear Rate =", min(-dVzdr), "[1/s]\n")

print("Max Viscosity =", max(eta), "[Pa*s];" , "Min Viscosity =", min(eta), "[Pa*s]\n")


eta_PL = K*(V_average/(2*R))**(n-1)*8**(n-1)*((3*n+1)/(4*n))**n
Re_PL = rho*V_average*2*R/eta_PL
L_e_pipe = Re_PL*2*R*0.06*1e6
print("Re_PL =", Re_PL, ";" , "Entrance Length =", L_e_pipe, "[micro-m]\n")

# validation of results
#print(dVzdr)
#print(tau_rz)
#print(eta)

#############################################################################################

# PANDAS DATAFRAME SECTION #

df = pd.read_csv(f'{data_file_name}.csv').dropna()

#sqrt(shearStress_XX^2+shearStress_YY^2+shearStress_ZZ^2+shearStress_XY^2+shearStress_YZ^2+shearStress_XZ^2)*rho
plot_cols_shearStress = [col for col in df.columns if 'shearStress' in col]
df_shearStress = (df[plot_cols_shearStress]*rho)**2
df_shearStress = np.array(df_shearStress.sum(axis=1)**(1/2))
xx_shearStress = np.linspace(1e-6, R*1e6, len(df_shearStress)) # convert unit to mircometer

plot_cols_shearRate = [col for col in df.columns if 'strainRate' in col]
df_shearRate = np.array((df[plot_cols_shearRate]))
xx_shearRate = np.linspace(1e-6, R*1e6, len(df_shearRate))

plot_cols_nu = [col for col in df.columns if 'nu' in col]
df_nu = np.array((df[plot_cols_nu]*rho))
xx_nu = np.linspace(1e-6, R*1e6, len(df_nu))

###########################
# power law fitting
from scipy.optimize import curve_fit
def power_law(gamma_dot_n, K, n):
    return K * np.power(gamma_dot_n, -n) # power law model

df_nu_fit = df_nu.ravel()
df_shearRate_fit = df_shearRate.ravel()
popt, pcov = curve_fit(power_law, df_shearRate_fit, df_nu_fit,) # min, max fiting y boundary considered

#residuals = df_nu_fit - power_law(df_shearRate_fit, *popt)
#ss_res = np.sum(residuals**2)
#ss_tot = np.sum((df_nu_fit-np.mean(df_nu_fit))**2)
#r_squared_scipy = 1 - (ss_res / ss_tot) # not vaild for non-linear fitting
###########################

plot_cols_U = [col for col in df.columns if 'U:' in col] # put : to aviod selecting grad(U)
df_U = (df[plot_cols_U])**2 # m/s
df_U = df_U.sum(axis=1)**(1/2)
xx_U = np.linspace(1e-6, R*1e6, len(df_U))

#############################################################################################
if plot_graphs == True:
    # PLOTTING SECTION #

    # plot settings
    linewidth = 2
    figsize_single = (6,6) # for signle plot use only
    figsize_double = (12,6) # for double plots use only
    title_font_weight = "bold"
    legend_size = {'size': 12}
    x_meter = x # x in meter in Vr(x_meter)
    x = np.linspace(1e-6, R*1e6, N) # covert meter to micrometer
    radius_name = "Variable Radius, $r$" #micrometer
    radius_name_w_unit = "Variable Radius, $r$ [$\\mu$m]"
    plot_title = True
    show_folder_images = False

    #############################################################################################
    # Velocity Profile
    fig1_name = "Velocity Profile, $V_z(r)$" # unit = [micrometer/s]
    fig1 = plt.figure(1, figsize=figsize_single)
    x_neg_R_pos_R = np.concatenate([-x[::-1],x[1:]])

    ###########################
    # Analytical Solution #
    plt.plot(x_neg_R_pos_R, np.concatenate([Vz(x_meter)[::-1]*1e3,Vz(x_meter)[1:]*1e3]), 'r', lw=linewidth, ls='--') # (e3)=millimeter/s, (e6)=micrometer/s
    # OpenFOAM Simulation #
    plt.plot(np.concatenate([-xx_U[::-1],xx_U[1:]]), np.concatenate([df_U[::-1]*1e3,df_U[1:]*1e3]), 'g', lw=linewidth, ls='-', alpha=0.7) # (e3)=millimeter/s, (e6)=micrometer/s
    ###########################

    #plt.autoscale(tight=True)
    #ax.set_xlim(0,R)
    #ax.set_ylim(bottom=0.)

    fig1.supxlabel(f"{radius_name_w_unit}")
    #fig1.supylabel(f"{fig1_name} [$\\mu$m/s]") # micrometer/s
    fig1.supylabel(f"{fig1_name} [mm/s]") # millimeter/s
    if plot_title == True:
        fig1.suptitle(f"{fig1_name} vs. {radius_name}", weight=title_font_weight)
    fig1.legend(['Analytical Solution', 'OpenFOAM Simulation'], loc='lower center', bbox_to_anchor=(0.5, 0.11), prop=legend_size)
    
    if save_graphs == True:
        plt.savefig('matplotlib_needle_figures/velocity_profile.png', dpi=dpi_save)

    #############################################################################################
    # Shear Rate
    fig2_name = "Shear Rate, $\\dot{\\gamma}$"
    fig2_name2 = "Shear Rate, $\\dot{\\gamma}$ [1/s]"
    fig2 = plt.figure(2, figsize=figsize_single)

    ###########################
    # Analytical Solution #
    plt.plot(x, -dVzdr, 'r', lw=linewidth, ls='--')
    # OpenFOAM Simulation #
    plt.plot(xx_shearRate, df_shearRate, 'g', lw=linewidth, ls='-', alpha=0.7)
    ###########################

    fig2.supxlabel(f"{radius_name_w_unit}")
    fig2.supylabel(f"{fig2_name} [1/s]")
    if plot_title == True:
        fig2.suptitle(f"{fig2_name} vs. {radius_name}", weight=title_font_weight)
    fig2.legend(['Analytical Solution', 'OpenFOAM Simulation'], loc='upper center', bbox_to_anchor=(0.33, 0.88), prop=legend_size)
    #fig2.legend(['Analytical Solution', 'OpenFOAM Simulation'], loc='lower center', bbox_to_anchor=(0.5, 0.15), prop=legend_size)
    ax = plt.gca()
    ax.set_yscale("log") # enable log scale on the y-axis (Shear Rate)

    if save_graphs == True:
        plt.savefig('matplotlib_needle_figures/shear_rate.png', dpi=dpi_save)

    #############################################################################################
    # Shear Stress
    fig3, axes = plt.subplots(1, 2, figsize=figsize_double)

    fig3a_name = "Shear Stress, $\\tau_{rz}$"
    fig3a_name2 = "Shear Stress, $\\tau_{rz}$ [Pa]"

    ###########################
    # Analytical Solution #
    axes[0].plot(x, tau_rz, 'r', lw=linewidth, ls='--', label="Analytical Solution")
    # OpenFOAM Simulation #
    axes[0].plot(xx_shearStress, df_shearStress, 'g', lw=linewidth, ls='-', label="OpenFOAM Simulation", alpha=0.7)
    ###########################

    if plot_title == True:
        axes[0].set_title(f"{fig3a_name} vs. {radius_name}", weight=title_font_weight)
    axes[0].set_xlabel(f"{radius_name_w_unit}")
    axes[0].set_ylabel(f"{fig3a_name2}")
    axes[0].legend(loc='upper left', prop=legend_size)

    fig3b_name = "Shear Stress, $\\tau_{rz}$ vs. Shear Rate, $\\dot{\\gamma}$"

    ###########################
    # Analytical Solution #
    axes[1].plot(-dVzdr, tau_rz, 'r', lw=linewidth, ls='--', label="Analytical Solution")
    # OpenFOAM Simulation #
    axes[1].plot(df_shearRate, df_shearStress, 'g', lw=linewidth, ls='-', label="OpenFOAM Simulation", alpha=0.7)
    ###########################

    #axes[1].ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    #axes[1].ticklabel_format(useMathText=True)
    if plot_title == True:
        axes[1].set_title(f"{fig3b_name}", weight="bold")
    axes[1].set_xlabel(f"{fig2_name2}")
    axes[1].set_ylabel(f"{fig3a_name2}")
    axes[1].legend(loc='upper left', prop=legend_size)

    if save_graphs == True:
        plt.savefig('matplotlib_needle_figures/shear_stress.png', dpi=dpi_save)

    #############################################################################################
    # Apparent Viscosity
    fig5, axes = plt.subplots(1, 2, figsize=figsize_double)

    fig5a_name = "Apparent Viscosity, $\\eta$"
    fig5a_name2 = "Apparent Viscosity, $\\eta$ [Pa$\\cdot{}$s]"

    ###########################
    # Analytical Solution #
    axes[0].plot(x, eta, 'r', lw=linewidth, ls='--', label="Analytical Solution", alpha=0.8)
    # OpenFOAM Simulation #
    axes[0].plot(xx_nu, df_nu, 'g', lw=linewidth, ls='-', label="OpenFOAM Simulation", alpha=0.7)
    ###########################

    if plot_title == True:
        axes[0].set_title(f"{fig5a_name} vs. {radius_name}", weight=title_font_weight)
    axes[0].set_xlabel(f"{radius_name_w_unit}")
    axes[0].set_ylabel(f"{fig5a_name2}")
    axes[0].legend(loc='upper right', prop=legend_size)

    fig5b_name = "Apparent Viscosity, $\\eta$ vs. Shear Rate, $\\dot{\\gamma}$"
    ###########################
    # Power Law Fitting Curve #
    power_law_df_shearRate_fit_popt = power_law(df_shearRate_fit, *popt)
    eta_equ_fit_power = f"$\\eta$ = {K}$\\dot\\gamma^{{{-n}}}$"
    axes[1].plot(-dVzdr, eta, 'r', lw=linewidth, ls='--', label=f"Power Law Fitting Curve, {eta_equ_fit_power}")
    #axes[1].set_ylim([0.2, 7000])

    # OpenFOAM Simulation #
    axes[1].plot(df_shearRate, df_nu, 'g', lw=2.5, ls='-', label="OpenFOAM Simulation", alpha=0.6)

    # Simulation Fitting Curve #
    eta_fit_name = "$\\eta_{fit}$"
    eta_equ_fit_scipy = eta_fit_name + " = {}$\\dot\\gamma^{{{}}}$".format(round(popt[0],3), round(popt[1]-1,3))
    axes[1].plot(df_shearRate_fit, power_law_df_shearRate_fit_popt, 'blue', lw=linewidth, ls='--', alpha=0.7,
    label=f"Simulation Fitting Curve, {eta_equ_fit_scipy}")
    ###########################
    # small plot inside
    # adjust accordingly for different data
    ax5s = fig5.add_axes([0.66, 0.41, 0.21, 0.21])
    plt.plot(-dVzdr, eta, 'r', lw=linewidth, ls='--', alpha=0.8)
    plt.plot(df_shearRate, df_nu, 'g', lw=2.5, ls='-', alpha=0.6)
    plt.plot(df_shearRate_fit, power_law(df_shearRate_fit, *popt), 'blue', lw=linewidth, ls='--', alpha=0.7)
    ax5s.set_ylim(0, 20)
    #ax5s.set_xlim(min(df_shearRate), max(df_shearRate))
    ax5s.set_xlabel("$\\dot\\gamma^{{{}}}$ [1/s]")
    plt.locator_params(axis='x', nbins=6)
    plt.locator_params(axis='y', nbins=7)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax5s.set_ylabel("$\\eta$ [Pa]")
    ax5s.ticklabel_format(useMathText=True)

    # Create a Rectangle patch

    rect = patches.Rectangle((df_shearRate_fit.min(), min(power_law(df_shearRate_fit, *popt))), 
        df_shearRate_fit.max()-df_shearRate_fit.min(), max(power_law_df_shearRate_fit_popt)-min(power_law_df_shearRate_fit_popt), 
        linewidth=1, edgecolor='k', facecolor='none')
    # Add the patch to the Axes
    axes[1].add_patch(rect)

    ###########################

    # Experimental Data #

    # NONE

    ###########################
    if plot_title == True:
        axes[1].set_title("{}".format(fig5b_name), weight="bold")
    axes[1].set_xlabel("{} [1/s]".format(fig2_name))
    axes[1].set_ylabel("{}".format(fig5a_name2))
    #axes[1].set_ylim([0.06, 0.18])
    axes[1].set_xscale("log")
    #axes[1].set_yscale("log")
    axes[1].legend(loc='upper right', prop={'size': 11-1})

    if save_graphs == True:
        plt.savefig('matplotlib_needle_figures/apparent_viscosity.png', dpi=dpi_save)
    # text box
    #props = dict(boxstyle='square', facecolor='none', edgecolor='none', ls='-', lw=linewidth, alpha=1)
    #axes[1].text(1.50, 0.65, eta_equ_fit, transform=ax.transAxes, fontsize=12, verticalalignment='top', bbox=props)

    #############################################################################################

    # Contour Plots

    cp = plt.figure(6, figsize=(7.5,6))
    [X,Y] = np.meshgrid(tau_rz,tau_rz)
    Z = np.sqrt(X**2 + Y**2)

    """
    # 1/4 cut contour
    # Analytical Solution #
    cpp = plt.contourf(X,Y,Z, levels=np.linspace(0,max(tau_rz),75), cmap=cm.coolwarm, antialiased=True, extend='neither')
    cbar = plt.colorbar(ticks=np.linspace(0,max(tau_rz),10))
    cp.suptitle("Contour Plot of {} within {}".format(fig3a_name, radius_name), weight=title_font_weight)
    cbar.set_label("{}".format(fig3a_name2))
    for c in cpp.collections:
        c.set_edgecolor("face")
    plt.axis('off')
    plt.savefig('matplotlib_needle_figures/1_4_cut_contour.png', dpi=dpi_save)
    """

    # full cut contour
    # Analytical Solution #
    tau_rz = np.array(tau_rz)
    tau_rz_rev = tau_rz[::-1]
    tau_rz = np.concatenate([-tau_rz_rev,tau_rz[1:]])
    [X,Y] = np.meshgrid(tau_rz,tau_rz)
    Z = np.sqrt(X**2 + Y**2)

    cpp = plt.contourf(x_neg_R_pos_R,x_neg_R_pos_R,Z, levels=np.linspace(0,max(tau_rz),75), cmap=cm.coolwarm, antialiased=True, extend='neither') # cmap=cm.rainbow, cmap=cm.jet
    cbar = plt.colorbar(ticks=np.linspace(0,max(tau_rz),10))
    cbar.ax.tick_params(labelsize=15)
    cbar.set_label(f"{fig3a_name2}",size=15)
    cp.supxlabel(f"{radius_name_w_unit}                  ") # centering xlabel
    cp.supylabel(f"{radius_name_w_unit}")
    #plt.axis('off')
    #plt.tight_layout()
    if plot_title == True:
        cp.suptitle(f"Contour Plot of {fig3a_name} (Analytical Solution) within Variable Radius, $-r$ to $r$", weight=title_font_weight, size=11)
    for c in cpp.collections:
        c.set_edgecolor("face")

    if save_graphs == True:
        plt.savefig('matplotlib_needle_figures/full_cut_contour.png', dpi=dpi_save)

    """
    # 3D
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    surf1 = ax.plot_surface(X, Y, Z, cmap=plt.cm.coolwarm, linewidth=0, antialiased=False)
    ax.contour(X, Y, Z, zdir='z', offset=-0.5)
    """

    #############################################################################################
    if show_folder_images == True:
        import matplotlib.image as mpimg
        plt.figure(7, figsize=(7.5,6)) # display paraview result to compare with the contour plot
        img_paraview = mpimg.imread('matplotlib_needle_figures/cross_section_paraview.png')
        imgplot = plt.imshow(img_paraview)
        plt.axis('off')
        plt.suptitle("Paraview Shear Stress Data, $y$ = 0.03 m", weight="bold")

        plt.figure(8, figsize=(7.5,6)) # residual plots
        img_paraview = mpimg.imread('matplotlib_needle_figures/residuals.png')
        imgplot = plt.imshow(img_paraview)
        plt.axis('off')
    #############################################################################################

    plt.show()
