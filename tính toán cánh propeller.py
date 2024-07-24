import math as m
#In this problem, all the number will show on SI system
#A propeller has a diameter of R (m) and rotates at a speed of rpm (rpm) when 
#the aircraft speed is V (m/s) at sea level in the standard atmosphere to produce 
#a thrust of N (N). Find 9 radial positions that divide the propeller’s swept 
#area into 8 equal parts when the hub - tip ratio is 0.2
#In this first case, I assume that the drag coefficient is very small
#density of the air at sea-level:
p0 = 1.225 #kg/m3
#Pressure of the air at sea-level:
P0 = 101325 #Pa (N/m2)
#temperature of the air at sea-level:
T0 = 288.15 #K
#we have the propeller diameter is D(m) so the raidus is equal to the half of 
#the propeller diameter 
D = float(input("input the diameter of the propeller cover (m), D= "))
R = round((D/2), 4)
T = float(input("input the thrust force created by the engine (N), T= "))
V = float(input("input the velocity of the aircraft at sea-level (m/s), V= "))
rpm = int(input("input the rpm (round per minute) of the flow aircraft "))
#hub ratio
hubratio = float(input("input the hub-tip ratio of the propeller "))
#Rh is the radius on the scale of the hub-tip ratio of the aircraft propeller
Rh = round((hubratio*R), 4)
#in this case, I will seperate the raidus onto 9 radial positions, so it will contains 
#8 difference radius of the propeller
#the first position is the mean radius also mean the R5
Rmean =round(m.sqrt((Rh**2+R**2)/2), 4)
R2 = round(m.sqrt((Rmean**2+Rh**2)/2), 4)
R1 = round(m.sqrt((R2**2+Rh**2)/2), 4)
R3 = round(m.sqrt((Rmean**2+R2**2)/2), 4)
R6 = round(m.sqrt((Rmean**2+R**2)/2), 4)
R5 = round(m.sqrt((R6**2+Rmean**2)/2), 4)
R7 = round(m.sqrt((R**2+R6**2)/2), 4)
omega = 2*m.pi*rpm/60
#priting all the radius position on display screen
#hide
print("The Rh radius position has the result of (m)", Rh)
print("The R1 radius position has the result of (m)", R1)
print("The R2 radius position has the result of (m)", R2)
print("The R3 radius position has the result of (m)", R3)
print("The Rmean radius position has the result of (m)", Rmean)
print("The R5 radius position has the result of (m)", R5)
print("The R6 radius position has the result of (m)", R6)
print("The R7 radius position has the result of (m)", R7)
print("The Rtip radius position has the result of (m)", R)
#We have to to use the BEMT (Blade element and Momentum theory) to find the axial inflow factor
#a and the angular factor inflow b
# a^2 + a - (T/2pV^2piR^2) = 0
# ax^2+ bx+ c = 0
#so to solve this equation we have to use the delta
#delta = b^2 - 4ac
#the pressure distribution and all the factor apply along with the propeller from the hub to the
#tip so that in this situation I want to find a average point at the center (mean radius) to
#demonstrate the axial inflow and angular inflow factor
#at first I will find the axial inflow factor a:
a1 = 1
b1 = 1
c1 = -T/(2*p0*(V**2)*m.pi*(R**2))
delta = b1**2 - 4*a1*c1
x1 = round((-b1+m.sqrt(delta))/(2*a1), 4)
x2 = round((-b1-m.sqrt(delta))/(2*a1), 4)
if delta > 1:
    xa = x1
else:
    xa = x2
print("axial inflow factor a is ", xa)
#next then I will find the angular inflow factor b:
#b^2  - b   +  T/(2pA(2*pi*rpm/60*R^2))=0
#ax^2 + bx  +  c
#the area here is the area of the propeller cover:
a2 = 1
b2 = -1
c2 = T/(2*p0*m.pi*R**2*(2*m.pi*rpm/60*Rmean)**2)
delta2 = b2**2 - 4*a2*c2
x3 = round((-b2+m.sqrt(delta2))/(2*a2), 4)
x4 = round((-b2-m.sqrt(delta2))/(2*a2), 4)
#the angular inflow factor b has 2 difference results but I must to choose which answer that
#will suit the condition that it can recall the same result of axial inflow factor a
test = x3**2 - x3 + xa*(1+xa)*V**2/(2*m.pi*rpm*R)**2
if 0.001 <= test <= 0.01:
    xb = x3
else:
    xb = x4
print("angular inflow factor b is ", xb)
#after finding the axial and angular inflow factor, I have to find the angle of the blade
#at the radius respectively
#the equation to find will base on the BEMT:
#tan(phi)=V*(1+a)/(omega*r*(1-b))
#I have to define the equation to find the tan(phi) for each radius
def tanphieq(r):
    return m.atan(V*(1+xa)/(omega*r*(1-xb)))/m.pi*  180
r_value = [Rh, R1, R2, R3, Rmean, R5, R6, R7, R]
for r in r_value:
    phi = tanphieq(r)
    print(f"for r = {r}, phi = {phi}")
#next the I have to find the local efficiency in the radial position of the blade
#the efficiency is equal to V*dF/(omega*dQ)
#the local efficiency is:
#dF = 4a(1+a)pV^2pirdr
#dQ = 4b(1+a)pVomegapir^3dr
def efficiency(r):
    return (xa/xb)*(V/omega*r)**2
r_value = [Rh, R1, R2, R3, Rmean, R5, R6, R7, R]
for r in r_value:
    eff = efficiency(r)
    print(f"for r = {r}, local efficiency = {eff}")

