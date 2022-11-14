import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mp
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

###############################################################################
# Single Pandas series OpenFOAM residuals plotting tool, KORIN, version 2.0.0 #
###############################################################################

# include residuals function in controlDict before simulation, OpenFOAM 9
"""
functions
{
	#includeFunc shearStress
	#includeFunc wallShearStress
	#includeFunc strainRate
	#includeFunc grad(U)
	
	residuals
    {
        type            		residuals;
        functionObjectLibs 		("libutilityFunctionObjects.so");
        enabled         		true;
        writeControl   		timeStep;
        writeInterval  		1;

        fields
        (
            U
            p
        );
    }
}
"""

file_name = "residuals"
series_name = "# Residuals         "

df = pd.read_csv(f'postProcessing/residuals/0/{file_name}.dat').dropna()
df = df[f'{series_name}'].str.split(pat="\t", expand=True)  #split the series into dataframe; select appropriate delimiter
df.columns = ['Time', 'Ux', 'Uy', 'Uz', 'p']                    #decide appropriate column names
df.drop(0,axis=0,inplace=True)                                  #cut first line which contains strings
df.to_csv(f'sorted_{file_name}.csv',index=False)                    #rewrite the dataframe into a new .csv file
df = pd.read_csv(f'sorted_{file_name}.csv')

df_logs = pd.read_csv('logs/contGlobal_0', sep='\t', names=['Iteration', 'Flux Field']).dropna()

print(df_logs)
print(df)

#############################################################################################
# Local Imbalance
fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(6,6))

linewidth = 1.15
#df = df[df.columns[1:]]
#df.plot(logx=True, xlim=(1,df.index.max()),style=['--','-','--','-'],alpha=0.8,
#   title='Local Imbalance',xlabel='Iteration Step',ylabel='Residual',ax=axes[0])

Ux = axes.plot(df['Time'], df['Ux'], ls='--', lw=linewidth, label='$U_x$', alpha=1)
Uy = axes.plot(df['Time'], df['Uy'], ls='-', lw=linewidth, label='$U_y$', alpha=0.5)
Uz = axes.plot(df['Time'], df['Uz'], ls='-', lw=linewidth, label='$U_z$', alpha=0.75)
p  = axes.plot(df['Time'], df['p'], ls='-', lw=linewidth, label='$p$', alpha=0.75)
axes.set_yscale("log")
axes.legend(loc='upper right', prop={'size': 12})
axes.set_xlim([1, len(df["Time"])])

axes.set_xlabel("Iteration Step")
axes.set_ylabel("Residual")
axes.set_title('Local Imbalance', weight='bold')

"""
contGlobal_0 = axes[1].plot(df_logs['Iteration'], df_logs['Flux Field'], 'm', label='Flux Field', alpha=0.75)
y000 = axes[1].axhline(y=0, color="#000000", ls='--', label='$y=0$', alpha=0.75)
#axes[1].set_xscale("log")
axes[1].legend(loc='upper right', prop={'size': 12})
axes[1].set_xlim([1, len(df_logs["Iteration"])])
axes[1].set_xlabel("Iteration Step")
axes[1].set_ylabel("Percentage Imbalance")
axes[1].set_title('Global Imbalance', weight='bold')
#df_logs.plot(logx=True, xlim=(1,df_logs.index.max()), color='m', x='Iteration', y='Flux Field', title='Global Imbalance', 
#   xlabel='Iteration Step',ylabel='Percentage Imbalance', ax=axes[1])
"""
#plt.savefig('residuals.png', dpi=600)

#############################################################################################
# Global Imbalance
#OpenFOAM --> foamLog foam.log
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, gridspec_kw={'width_ratios': (1, 2)}, figsize=(6,6))

###########################
ax1.set_xlim(0, 60)  # outliers only
ax2.set_xlim(600, len(df_logs["Iteration"]))  # most of the data
###########################

fig.subplots_adjust(wspace=0.15)
ax1.plot(df_logs['Iteration'], df_logs['Flux Field'], 'm', lw=linewidth+0.25, label='Flux Field', alpha=0.75)
ax2.plot(df_logs['Iteration'], df_logs['Flux Field'], 'm', lw=linewidth+0.25, label='Flux Field', alpha=0.75)
ax1.axhline(y=0, color="#000000", ls='--', lw=linewidth+0.25, label='$y=0$', alpha=0.75)
ax2.axhline(y=0, color="#000000", ls='--', lw=linewidth+0.25, label='$y=0$', alpha=0.75)
ax2.legend(loc='upper right', prop={'size': 12})
ax1.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax2.tick_params(axis='y', length=0)
ax1.set_ylabel("Percentage Imbalance")
ax2.set_xlabel("Iteration Step                                         ")
ax2.set_title('Global Imbalance                                ', weight='bold')

d = 1  # proportion of vertical to horizontal extent of the slanted line
kwargs = dict(marker=[(-1, -d), (1, d)], markersize=8,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
ax1.plot([1, 1], [1, 0], transform=ax1.transAxes, **kwargs)
ax2.plot([0, 0], [1, 0], transform=ax2.transAxes, **kwargs)

#plt.savefig('global_imbalance.png', dpi=600)

#############################################################################################

plt.show()
