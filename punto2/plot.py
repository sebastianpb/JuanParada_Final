from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np

data=np.loadtxt( 'datos_observacionales.dat' )
t,x,y,z=data.T
sigma=y.std()

def differential_eqs(w, t, p):
    """
    Sistema de ecuaciones diferenciales acopladas
        w :  [x,y,z]
        t :  tiempo
        p :  vector de parámetros p = [sigma,rho,beta]
    """
    x,y,z = w
    sigma,rho,beta = p

    f = [sigma*(y-x),
         (x*(rho-z)) - y,
         (x*y) - (beta*z)]
    return f

    
'''
p=[1,1,1]
print(w0)
sols = odeint(differential_eqs,w0,t,args=(p,))
    
print(sols)
plt.plot(t,sols)
'''
def model(t,w0,params):
    sols=odeint(differential_eqs,w0,t,args=(params,))
    return np.array(sols[:,1])  # Retorna solo y

def log_likelihood(t,y,w0,sigma,params):
    return -0.5*np.sum(((y-model(t,w0,params))/(sigma))**2)


def gradient_log_pdf_to_sample(x,y, w0,sigma,params):
    dx=1e-5
    grad=np.zeros(len(params))
    for i in range(len(params)):
        delta=np.zeros(len(params))
        delta[i]=dx
        f1=log_likelihood(x,y,w0,sigma,params+delta)
        f2=log_likelihood(x,y,w0,sigma,params)
        grad[i]=(f1-f2)/dx
    return(grad)

def leapfrog(x,y, w0,sigma, params ,p, delta_t=1E-1, niter=5):
    '''
    Integracion tipo leapfrog. 
        `params` representa las posiciones (i.e. los parametros).
        `p` representa el momentum asociado a los parametros.
    '''
    m = 100.0
    q_new = params.copy()
    p_new = p.copy()
    for i in range(niter):
        p_new = p_new + 0.5 * delta_t * gradient_log_pdf_to_sample(x,y,w0,sigma,q_new) #kick
        q_new = q_new + delta_t * (p_new/m) #drift
        p_new = p_new + 0.5 * delta_t * gradient_log_pdf_to_sample(x,y,w0,sigma,q_new) #kick
    return q_new, p_new

    
def H(x,y,w0,sigma, params,p):
    m=100
    K = 0.5 * (np.sum(p**2))/m
    U = -log_likelihood(x,y,w0,sigma,params)
    return K + U

def MCMC(t,y,w0,sigma,nsteps=1000):
    n_param=3 # Numero de parametros!!
    q = np.zeros((nsteps,n_param))
    p = np.zeros((nsteps,n_param))
    p[0,:] = np.random.normal(0,1,n_param)
    q[0,:] = np.random.normal(0,1,n_param)
    #sigma = 0.1
    for i in range(1,nsteps):

        q_new, p_new = leapfrog(t,y,w0,sigma,q[i-1,:],p[i-1]) # la propuesta se hace con leapfrog
        p_new = -p_new #negamos a p para que la propuesta sea simetrica.
        
        E_new = H(x,y,w0,sigma,q_new,p_new) # En lugar de evaluar la pdf se evalua la energia.
        E_old = H(x,y,w0,sigma,q[i-1,:],p[i-1])
        alpha = min(1.0,np.exp(-(E_new - E_old))) # Se comparan las dos energias
        
        beta = np.random.random()
        if beta < alpha:
            q[i,:] = q_new
        else:
            q[i,:] = q[i-1]
    return q

# plt.errorbar(t, y, yerr=sigma, fmt='o')
# initial condition
w0 = list(data[0,1:])

param_chain=MCMC(t,y,w0,sigma,500)
n_param  = len(param_chain[0])
best = []
for i in range(n_param):
    best.append(np.mean(param_chain[:,i]))
	

def all_model(t,w0,params):
    sols=odeint(differential_eqs,w0,t,args=(params,))
    return np.array(sols)  # Retorna solo y


t_model = np.linspace(min(t), max(t), 100)
sols_model = all_model(t_model, w0,best)
x_model,y_model,z_model=sols_model.T

plt.errorbar(t,x, yerr=sigma, fmt='o', label='x_obs')
plt.errorbar(t,y, yerr=sigma, fmt='o', label='y_obs')
plt.errorbar(t,z, yerr=sigma, fmt='o', label='z_obs')
plt.plot(t_model, x_model, label='x_model')
plt.plot(t_model, y_model, label='y_model')
plt.plot(t_model, z_model, label='z_model')
plt.xlabel("t")
plt.legend()
plt.savefig("f1.pdf")


chain_sigma=param_chain[:,0]
chain_rho=param_chain[:,1]
chain_beta=param_chain[:,2]

def title(chain):
    v_medio  = np.mean(chain)
    v_sigma = np.std(chain)
    title = "$v = {:.5f} \pm {:.5f} $ km/s".format(v_medio, v_sigma)

    return title

fig, axs = plt.subplots(1, 3, figsize=(10,5))
axs[0].hist(chain_sigma, bins=50, density=True)
axs[0].set_title(title(chain_sigma))
axs[1].hist(chain_rho, bins=50, density=True)
axs[1].set_title(title(chain_rho))
axs[2].hist(chain_beta, bins=50, density=True)
axs[2].set_title(title(chain_beta))
plt.savefig("f2.pdf")