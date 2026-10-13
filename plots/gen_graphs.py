import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

N = [400, 600, 800, 1000, 1200]

neopt = [0.644, 2.173, 5.649, 9.606, 18.399]
blas  = [0.00860, 0.01792, 0.03843, 0.07123, 0.12376]
opt   = [0.209, 0.734, 1.694, 3.267, 5.554]

# --- Graf 1: Timp de executie (scara liniara) ---
fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(N, neopt, 'o-', color='tab:red',   label='neopt',  linewidth=2, markersize=6)
ax.plot(N, opt,   's-', color='tab:orange', label='opt_m',  linewidth=2, markersize=6)
ax.plot(N, blas,  '^-', color='tab:blue',  label='blas',   linewidth=2, markersize=6)
ax.set_xlabel('Dimensiunea matricei N', fontsize=12)
ax.set_ylabel('Timp de executie (secunde)', fontsize=12)
ax.set_title('Timp de executie vs. N (partitia Haswell, -O0)', fontsize=13)
ax.legend(fontsize=11)
ax.grid(True, linestyle='--', alpha=0.5)
ax.set_xticks(N)
plt.tight_layout()
plt.savefig('timp_liniar.png', dpi=150)
plt.close()
print('Saved timp_liniar.png')

# --- Graf 2: Timp de executie (scara logaritmica) ---
fig, ax = plt.subplots(figsize=(9, 5))
ax.semilogy(N, neopt, 'o-', color='tab:red',   label='neopt',  linewidth=2, markersize=6)
ax.semilogy(N, opt,   's-', color='tab:orange', label='opt_m',  linewidth=2, markersize=6)
ax.semilogy(N, blas,  '^-', color='tab:blue',  label='blas',   linewidth=2, markersize=6)
ax.set_xlabel('Dimensiunea matricei N', fontsize=12)
ax.set_ylabel('Timp de executie (secunde, scala log)', fontsize=12)
ax.set_title('Timp de executie vs. N — scala logaritmica (partitia Haswell, -O0)', fontsize=13)
ax.legend(fontsize=11)
ax.grid(True, which='both', linestyle='--', alpha=0.5)
ax.set_xticks(N)
plt.tight_layout()
plt.savefig('timp_log.png', dpi=150)
plt.close()
print('Saved timp_log.png')

# --- Graf 3: Speedup fata de neopt ---
speedup_blas = [neopt[i] / blas[i] for i in range(len(N))]
speedup_opt  = [neopt[i] / opt[i]  for i in range(len(N))]

fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(N, speedup_blas, '^-', color='tab:blue',   label='speedup blas vs neopt',  linewidth=2, markersize=6)
ax.plot(N, speedup_opt,  's-', color='tab:orange', label='speedup opt_m vs neopt', linewidth=2, markersize=6)
ax.set_xlabel('Dimensiunea matricei N', fontsize=12)
ax.set_ylabel('Speedup (x ori mai rapid decat neopt)', fontsize=12)
ax.set_title('Speedup relativ fata de neopt (partitia Haswell, -O0)', fontsize=13)
ax.legend(fontsize=11)
ax.grid(True, linestyle='--', alpha=0.5)
ax.set_xticks(N)
for i, (xv, yb, yo) in enumerate(zip(N, speedup_blas, speedup_opt)):
    ax.annotate(f'{yb:.0f}x', (xv, yb), textcoords='offset points',
                xytext=(0, 6), ha='center', fontsize=9, color='tab:blue')
    ax.annotate(f'{yo:.1f}x', (xv, yo), textcoords='offset points',
                xytext=(0, 6), ha='center', fontsize=9, color='tab:orange')
plt.tight_layout()
plt.savefig('speedup.png', dpi=150)
plt.close()
print('Saved speedup.png')

# --- Graf 4: Cachegrind - D1 miss rate ---
labels = ['neopt', 'opt_m', 'blas']
d1_miss  = [4.4, 1.3, 5.3]   # D1 miss rate % din fisierele .cache (neopt, opt_m, blas)
lld_miss = [0.0, 0.0, 0.4]

x = np.arange(len(labels))
width = 0.35
fig, ax = plt.subplots(figsize=(7, 5))
bars1 = ax.bar(x - width/2, d1_miss, width, label='D1 miss rate (%)', color='tab:red', alpha=0.8)
bars2 = ax.bar(x + width/2, lld_miss, width, label='LLd miss rate (%)', color='tab:purple', alpha=0.8)
ax.set_xlabel('Implementare', fontsize=12)
ax.set_ylabel('Miss rate (%)', fontsize=12)
ax.set_title('Cache miss rate (Cachegrind, N=400, Haswell)', fontsize=13)
ax.set_xticks(x)
ax.set_xticklabels(labels, fontsize=11)
ax.legend(fontsize=11)
ax.grid(True, axis='y', linestyle='--', alpha=0.5)
for bar in bars1:
    h = bar.get_height()
    if h > 0:
        ax.annotate(f'{h}%', xy=(bar.get_x() + bar.get_width()/2, h),
                    xytext=(0, 3), textcoords='offset points', ha='center', fontsize=9)
for bar in bars2:
    h = bar.get_height()
    if h > 0:
        ax.annotate(f'{h}%', xy=(bar.get_x() + bar.get_width()/2, h),
                    xytext=(0, 3), textcoords='offset points', ha='center', fontsize=9)
plt.tight_layout()
plt.savefig('cache_miss.png', dpi=150)
plt.close()
print('Saved cache_miss.png')
