"""Plot the SpaceX Roadster in the Heliocentric frame

The ephemerides for the Roadster can be found at the following locations:

-143205 : JPL Horizons

https://www.projectpluto.com/temp/spacex.htm
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from astro import planets, constants, planets, time, kepler

au2km = constants.au2km
deg2rad = constants.deg2rad
sec2day = constants.sec2day
day2sec = constants.day2sec

# roadster heliocentric orbital elements
def roadster_coe(jd):
    """Output the COES of the roadster at the given JD
    """
    epoch = 2458190.500000000
    ecc = 2.560214125839113E-01  # eccentricy 
    per = 9.860631722474498E-01 * au2km # periapsis
    inc = 1.078548455272089E+00 * deg2rad # inclination
    raan = 3.171871226218810E+02 * deg2rad #  
    arg_p = 1.774087594994931E+02 * deg2rad 
    jd_per =  2458153.531206728425
    mean_motion = 6.459332187157216E-01 * deg2rad * sec2day 
    mean_anom = 2.387937163002724E+01 * deg2rad
    E, nu, count = kepler.kepler_eq_E(mean_anom, ecc)

    nu = 4.031991175737090E+01 * deg2rad
    sma = 1.325391871387247E+00 * au2km 
    apo = 1.664720570527044E+00 * au2km 
    period = 5.573331570030891E+02 * day2sec
    
    # compute some other elements
    p = kepler.semilatus_rectum(sma, ecc)
    deltat = (jd - epoch) * day2sec
    # propogate from epoch to JD and find nu
    Ef, Mf, nuf = kepler.tof_delta_t(p, ecc, constants.sun.mu, nu, deltat)

    roadster_coe = planets.COE(p=p, ecc=ecc, inc=inc, raan=raan,
                               argp=arg_p, nu=nuf)

    # compute the vector in J2000 Ecliptic frame
    r_ecliptic, v_ecliptic, r_pqw, v_pqw = kepler.coe2rv(p, ecc, inc,
                                                         raan, arg_p, nuf,
                                                         constants.sun.mu)

    return roadster_coe, r_ecliptic, v_ecliptic

# scale everything by AU2KM

jd_start = time.date2jd(2018, 3, 1, 0,0, 0)[0]
jd_end = time.date2jd(2028, 3, 1, 0, 0, 0)[0]

jd_span = np.arange(jd_start, jd_end, 1)

_, roadster_ecliptic, roadster_ecliptic = roadster_coe(jd_start)
mercury_coe, mercury_ecliptic, _, _ ,  _ = planets.planet_coe(jd_start, 0)
venus_coe, venus_ecliptic, _, _ ,  _ = planets.planet_coe(jd_start, 1)
earth_coe, earth_ecliptic, _, _, _= planets.planet_coe(jd_start, 2)
mars_coe, mars_ecliptic, _, _, _ = planets.planet_coe(jd_start, 3)

# initialize all the plot elements (orbit position and conic orbits)
fig, ax = plt.subplots()
ax.set_xlim(-5*au2km, 5*au2km)
ax.set_ylim(-5*au2km, 5*au2km)

rp, = ax.plot(roadster_ecliptic[0], roadster_ecliptic[1], marker='o', color='r')

mcp, = ax.plot(mercury_ecliptic[0], mercury_ecliptic[1], marker='o', color='c')
vp, = ax.plot(venus_ecliptic[0], venus_ecliptic[1], marker='o', color='m')
ep, = ax.plot(earth_ecliptic[0], earth_ecliptic[1],marker='o', color='g')
mp, = ax.plot(mars_ecliptic[0], mars_ecliptic[1],marker='x', color='b')

# draw conic orbits
mercury_inertial, _, _, _, _, _  = kepler.conic_orbit(mercury_coe.p, mercury_coe.ecc, mercury_coe.inc, mercury_coe.raan, mercury_coe.argp, 0, 2*np.pi, mu=constants.sun.mu)
venus_inertial, _, _, _, _, _  = kepler.conic_orbit(venus_coe.p, venus_coe.ecc, venus_coe.inc, venus_coe.raan, venus_coe.argp, 0, 2*np.pi, mu=constants.sun.mu)
earth_inertial, _, _, _, _, _  = kepler.conic_orbit(earth_coe.p, earth_coe.ecc, earth_coe.inc, earth_coe.raan, earth_coe.argp, 0, 2*np.pi, mu=constants.sun.mu)
mars_inertial, _, _, _, _, _  = kepler.conic_orbit(mars_coe.p, mars_coe.ecc, mars_coe.inc, mars_coe.raan, mars_coe.argp, 0, 2*np.pi, mu=constants.sun.mu)


ax.plot(mercury_inertial[:,0],mercury_inertial[:,1], color='c')
ax.plot(venus_inertial[:,0],venus_inertial[:,1], color='m')
ax.plot(earth_inertial[:,0],earth_inertial[:,1], color='green')
ax.plot(mars_inertial[:,0],mars_inertial[:,1], color='red')

ax.axis('equal')
def init():
    rp.set_data([], [])
    mcp.set_data([], [])
    vp.set_data([], [])
    ep.set_data([], [])
    mp.set_data([], [])

    # ax.plot(state_eci[:,0], state_eci[:,1], color='red')
    
    return rp, mcp, vp, ep, mp

def animate(jd):
    # compute the new position
    _, roadster_ecliptic, _ = roadster_coe(jd)
    _, mercury_ecliptic, _, _ ,  _ = planets.planet_coe(jd, 0)
    _, venus_ecliptic, _, _ ,  _ = planets.planet_coe(jd, 1)
    _, earth_ecliptic, _, _, _= planets.planet_coe(jd, 2)
    _, mars_ecliptic, _, _, _ = planets.planet_coe(jd, 3)

    rp.set_data(roadster_ecliptic[0], roadster_ecliptic[1])
    mcp.set_data(mercury_ecliptic[0], mercury_ecliptic[1])
    vp.set_data(venus_ecliptic[0], venus_ecliptic[1])
    ep.set_data(earth_ecliptic[0], earth_ecliptic[1])
    mp.set_data(mars_ecliptic[0], mars_ecliptic[1])
    
    return rp, mcp, vp, ep, mp

# get the elements for all the planets

# plot them all in an animation using matplotlib
ani = animation.FuncAnimation(fig, animate, jd_span, interval=10, blit=False, init_func=init)
# loop over julian date
plt.show()
