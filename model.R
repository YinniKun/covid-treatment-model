library(deSolve)

#parameter assignment - Data from COVID paper
# target cell model
Ro = 3.79 # 3.79
d = 0.001 # 0.001 rate of death of uninfected cells
k = 0.55 # 0.55 virus infect uninfected cells at rate k
delta = 0.11 # 0.11 infected cell death rate 
p = 0.24 # 0.24 rate of production of virus from infected cells
c = 5.36  # 5.36 death rate of virus
lambda = Ro*d*delta*c/(k*p)  # Ro*d*delta*c/(k*p); rate of production of uninfected cells

# target cell model with immune response
omiga = 4 #transition rate between I1 and I2, see influenza paper for this parameter
s = 0.1 #rate of secretion of IFN (randomly picked)
b = 0.1 #rate of degradation of IFN (randomly picked)
ek = 1 #effectiveness of blocking infection
ep = 1 #effectiveness of blocking production of virus
ew = 1 #effectiveness of IFN blocking transition 
tau = 1 #time delay before the INF production kicks in (randomly picked)

# drug X
mx = 1 #production
nx = 1 #elimination
Kx = 1 #equilibrium constant

#drug y
my = 1
ny = 1
Ky = 1

#drug Z
mz = 1
nz = 1
Kz = 1

#drug A
ma = 1
na = 1
ka =1 

#initial condition
Io = 0
Co = 25 - Io
Vo = 0.061


#######

## Original Model
model1 <- function (t, x, params) {
  C <- x[1] #target cell
  I <- x[2] #infected cell
  V <- x[3] #viral load
  
  dCdt = lambda - d*C - k*V*C 
  dIdt = k*V*C - delta * I
  dVdt = p*I-c*V
  
  dxdt <- c(dCdt, dIdt, dVdt)
  list(dxdt)
}

## Model with immune
derivs <- function(t, y, parms) {
  if (t < tau) {
    dCdt <- lambda - k*y[1]*y[4] - d*y[1] #uninfected cell
    dIadt <- k*y[1]*y[4] - omiga*y[2] #infected but not virus producing cells
    dIbdt <- omiga*y[2] - delta *y[3] #infected and virus producing cells
    dVdt <- p*y[3] - c*y[4] #viral load
    dNdt <- 0 #interferon produced by I2
    }
  else { #C=y[1],Ia=y[2],Ib=y[3],V=y[4],N=y[5]
    lag <- lagvalue(t - tau)
    
    dCdt = lambda - k*y[1]*y[4] - d*y[1]
    dIadt = k*y[1]*y[4] - (omiga/(1+ew*y[5]))*y[2]
    dIbdt = (omiga/(1+ew*y[5]))*y[2] - delta*y[2]
    dVdt = p/(1+ep*y[5])*y[3]-c*y[4]
    dNdt = s*lag[3] - b*y[5]
    }
  list(c(dCdt, dIadt, dIbdt, dVdt, dNdt))
}

## Extension of the model
## Model with Drug X, assume drug is taken after tau time units
extended <- function(t, y, parms) {
  if (t < tau) {
    dCdt <- lambda - k*y[1]*y[4] - d*y[1] #uninfected cell
    dIadt <- k*y[1]*y[4]-omiga*y[2] #infected but not virus producing cells
    dIbdt <- omiga*y[2] - delta *y[3] #infected and virus producing cells
    dVdt <- p*y[3] - c*y[4] #viral load
    dNdt <- 0 #interferon produced by I2
    dXdt <- 0
  }
  else { #C=y[1],Ia=y[2],Ib=y[3],V=y[4],N=y[5]
    lag <- lagvalue(t - tau)
    
    dCdt = lambda - k*y[1]*y[4] - d*y[1]
    dIadt = (k*y[1]*y[4])/(1+y[6]/Kx) - (omiga/(1+ew*y[5]))*y[2]
    dIbdt = (omiga/(1+ew*y[5]))*y[2] - delta*y[2]
    dVdt = p/(1+ep*y[5])*y[3]-c*y[4]
    dNdt = s*lag[3] - b*y[5]
    dXdt = mx - nx*y[6]
  }
  list(c(dCdt, dIadt, dIbdt, dVdt, dNdt, dXdt))
}

# Drug Y, assume drug is taken after tau time units
extended2 <- function(t, y, parms) {
  if (t < tau) {
    dCdt <- lambda - k*y[1]*y[4] - d*y[1] #uninfected cell
    dIadt <- k*y[1]*y[4] - omiga*y[2] #infected but not virus producing cells
    dIbdt <- omiga*y[2] - delta *y[3] #infected and virus producing cells
    dVdt <- p*y[3] - c*y[4] #viral load
    dNdt <- 0 #interferon produced by I2
    dYdt <- 0
  }
  else { #C=y[1],Ia=y[2],Ib=y[3],V=y[4],N=y[5]
    lag <- lagvalue(t - tau)

    dCdt = lambda - k*y[1]*y[4] - d*y[1]
    dIadt = k*y[1]*y[4] - (omiga/(1+ew*y[5]))*y[2]
    dIbdt = (omiga/(1+ew*y[5]))*y[2] - delta*y[2]
    dVdt = p/(1+ep*y[5])*y[3]-c*y[4]
    dNdt = s*lag[3] - b*y[5] + Ky*y[6]
    dYdt = my - ny*y[6]
  }
  list(c(dCdt, dIadt, dIbdt, dVdt, dNdt, dYdt))
}

# Drug Z, assume drug is taken after tau time units
extended3 <- function(t, y, parms) {
  if (t < tau) {
    dCdt <- lambda - k*y[1]*y[4] - d*y[1] #uninfected cell
    dIadt <- k*y[1]*y[4] - omiga*y[2] #infected but not virus producing cells
    dIbdt <- omiga*y[2] - delta *y[3] #infected and virus producing cells
    dVdt <- p*y[3] - c*y[4] #viral load
    dNdt <- 0 #interferon produced by I2
    dZdt <- 0
  }
  else { #C=y[1],Ia=y[2],Ib=y[3],V=y[4],N=y[5]
    lag <- lagvalue(t - tau)
    
    dCdt = lambda - k*y[1]*y[4] - d*y[1]
    dIadt = k*y[1]*y[4] - (omiga/(1+ew*y[5]))*y[2]
    dIbdt = (omiga/(1+ew*y[5]))*y[2] - delta*y[2]
    dVdt = p/(1+ep*y[5])*y[3]-c*y[4]
    dNdt = s*lag[3] - (b*y[5])/(1+(y[6]/Kz))
    dZdt = mz - nz*y[6]
  }
  list(c(dCdt, dIadt, dIbdt, dVdt, dNdt, dZdt))
}

# Drug A, assume drug is taken after tau time units
extended4 <- function(t, y, parms) {
  if (t < tau) {
    dCdt <- lambda - k*y[1]*y[4] - d*y[1] #uninfected cell
    dIadt <- k*y[1]*y[4] - omiga*y[2] #infected but not virus producing cells
    dIbdt <- omiga*y[2] - delta *y[3] #infected and virus producing cells
    dVdt <- p*y[3] - c*y[4] #viral load
    dNdt <- 0 #interferon produced by I2
    dAdt <- 0
  }
  else { #C=y[1],Ia=y[2],Ib=y[3],V=y[4],N=y[5]
    lag <- lagvalue(t - tau)
    
    dCdt = lambda - k*y[1]*y[4] - d*y[1]
    dIadt = k*y[1]*y[4] - (omiga/(1+ew*y[5]))*y[2] - ka*y[6]*y[2]
    dIbdt = (omiga/(1+ew*y[5]))*y[2] - delta*y[2]
    dVdt = p/(1+ep*y[5])*y[3]-c*y[4]
    dNdt = s*lag[3] - b*y[5]
    dAdt = ma - na*y[6] - ka*y[6]*y[2]
  }
  list(c(dCdt, dIadt, dIbdt, dVdt, dNdt, dAdt))
}

### Parameters
times <- seq(from=0, to=150, by=0.01)

x_initial <- c(C=Co, I=Io, V=Vo) #initial for target cell 
yinit <- c(Co,0,0,Vo,0) #initial for target cell with immune
yinit2 <- c(Co,0,0,Vo,0,0) #initial for drug x, y, and z

# generate the list
out <- as.data.frame(ode(func = model1, y=x_initial, parms = parms, times = times))
yout <- dede(y = yinit, times = times, func = derivs, parms = NULL)
yout2 <- dede(y = yinit2, times = times, func = extended, parms = NULL)
yout3 <- dede(y = yinit2, times = times, func = extended2, parms = NULL)
yout4 <- dede(y = yinit2, times = times, func = extended3, parms = NULL)
#yout5 <- dede(y = yinit2, times = times, func = extended4, parms = NULL)

# Plot for target cell
plot(out$C~out$time, type = "l", col = "red", xlab = "Time (in days)", ylab = "Concentration", main = "Simulation for Target Cell Model", ylim = c(0,0.7))
points(out$I~out$time, type = "l", col = "blue")
points(out$V~out$time, type = "l", col = "green")

# plot for target cell with immune response
# plot(yout[,2]~yout[,1], type = "l", col = "red", xlab = "Time (in day)", ylab = "Concentration", main = "Simulation for Target Cell Model with IIR", ylim = c(0,25))
# points(yout[,3]~yout[,1], type = "l", col = "blue")
# points(yout[,4]~yout[,1], type = "l", col = "green")
# points(yout[,5]~yout[,1], type = "l", col = "black")
# points(yout[,6]~yout[,1], type = "l", col = "orange")

# plot for drug x
# plot(yout2[,2]~yout2[,1], type = "l", col = "red", xlab = "Time (in day)", ylab = "Concentration", main = "Simulation for Target Cell Model with IIR and Drug X", ylim = c(0,25))
# points(yout2[,3]~yout2[,1], type = "l", col = "blue")
# points(yout2[,4]~yout2[,1], type = "l", col = "green")
# points(yout2[,5]~yout2[,1], type = "l", col = "black")
# points(yout2[,6]~yout2[,1], type = "l", col = "orange")
## points(yout2[,7]~yout2[,1], type = "l", col = "purple")

# plot for drug y
# plot(yout3[,2]~yout3[,1], type = "l", col = "red", xlab = "Time (in day)", ylab = "Concentration", main = "Simulation for Target Cell Model with IIR and Drug Y", ylim = c(0,25))
# points(yout3[,3]~yout3[,1], type = "l", col = "blue")
# points(yout3[,4]~yout3[,1], type = "l", col = "green")
# points(yout3[,5]~yout3[,1], type = "l", col = "black")
# points(yout3[,6]~yout3[,1], type = "l", col = "orange")
## points(yout3[,7]~yout3[,1], type = "l", col = "purple")

# plot for drug z
# plot(yout4[,2]~yout4[,1], type = "l", col = "red", xlab = "Time (in day)", ylab = "Concentration", main = "Simulation for Target Cell Model with IIR and Drug Z", ylim = c(0,25))
# points(yout4[,3]~yout4[,1], type = "l", col = "blue")
# points(yout4[,4]~yout4[,1], type = "l", col = "green")
# points(yout4[,5]~yout4[,1], type = "l", col = "black")
# points(yout4[,6]~yout4[,1], type = "l", col = "orange")
## points(yout4[,7]~yout4[,1], type = "l", col = "purple")

# plot for drug A
# plot(yout5[,2]~yout5[,1], type = "l", col = "red", xlab = "Time (in day)", ylab = "Concentration", main = "Simulation for Target Cell Model with IIR and Drug A", ylim = c(0,0.7))
# points(yout5[,3]~yout5[,1], type = "l", col = "blue")
# points(yout5[,4]~yout5[,1], type = "l", col = "green")
# points(yout5[,5]~yout5[,1], type = "l", col = "black")
# points(yout5[,6]~yout5[,1], type = "l", col = "orange")
##points(yout5[,7]~yout5[,1], type = "l", col = "purple")

# legends
legend(100,0.4,legend = c("Healthy Cells", "Infected Cells", "Viral Load"), col = c("red", "blue", "green"), pch=20)
#legend(100,10,legend = c("Healthy Cells", "Infected Cells 1", "Infected Cells 2", "Viral Load", "IFN"), col = c("red", "blue", "green", "black", "orange"), pch=20)
