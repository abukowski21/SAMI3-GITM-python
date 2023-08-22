import cartopy.crs as ccrs
from cartopy.feature.nightshade import Nightshade
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from utility_programs.utils import ut_to_lt
from matplotlib.dates import DateFormatter



def timedeltatotime(td, secs_back=False):
    td_in_seconds = td.total_seconds()
    hours, remainder = divmod(td_in_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)

    hours = int(hours)
    minutes = int(minutes)
    seconds = int(seconds)

    if hours <= -1:
        hours += 1
        if secs_back:
            tstring = '-%s:%s:%s' % (str(np.abs(hours)).rjust(2, '0'),
                                     str(minutes).rjust(2, '0'), str(seconds).rjust(2, '0'))
        else:
            tstring = '-%s:%s' % (str(np.abs(hours)).rjust(2,
                                  '0'), str(minutes).rjust(2, '0'))
    else:
        if secs_back:
            tstring = '%s:%s:%s' % (str(hours).rjust(2, '0'), str(
                minutes).rjust(2, '0'), str(seconds).rjust(2, '0'))
        else:
            tstring = '%s:%s' % (str(hours).rjust(
                2, '0'), str(minutes).rjust(2, '0'))

    return tstring


def fig1(ds,
         vlims,
         plot_start,
         plot_delta,
         storm_start,
         at_alt=None,
         numplots=6,
         suptitle_colname='Total Neutral Mass Density',
         mag_eq=None,
         ):

    fig, axs = plt.subplots(3, 2, figsize=(15, 10), subplot_kw={
        'projection': ccrs.PlateCarree()},
        sharex=True,
        sharey=True,)

    # plotting...
    # Setup some things...
    tnow = plot_start
    t = 0

    axes = axs.flatten()

    while tnow < plot_start + plot_delta * numplots:
        im = ds.sel(time=tnow, method='nearest').plot(x='lon', transform=ccrs.PlateCarree(),
                                                      vmin=-vlims,
                                                      vmax=vlims,
                                                      ax=axes[t],
                                                      add_colorbar=False,
                                                      cmap='bwr'
                                                      )

        axes[t].set_title("%s From Storm Onset" %
                          timedeltatotime(tnow-storm_start), fontsize='x-large')
        axes[t].coastlines(zorder=3, color='black', alpha=1)
        lines = axes[t].gridlines(color='black', linestyle='--',
                                  alpha=0.6, )
        lines.bottom_labels = True if t == 4 or t == 5 else False
        lines.left_labels = True if t == 0 or t == 2 or t == 4 else False

        axes[t].add_feature(Nightshade(tnow, alpha=0.21))
        
        if mag_eq is not None:
            axes[t].plot(mag_eq['glon'], mag_eq['mlat'], transform=ccrs.PlateCarree(),
                         linestyle='--', color='k')
        
        
        
        tnow += plot_delta
        t += 1

    fig.supxlabel('Longitude (Degrees East)', fontsize='x-large')
    fig.supylabel('Latitude (Degrees North)', fontsize='x-large')

    fig.suptitle("Perturbation Total Neutral Mass Density from GITM at %ikm" % at_alt,
                 fontsize='xx-large', fontweight='heavy')

    fig.tight_layout()

    cbar = fig.colorbar(im, ax=axs.ravel().tolist(), orientation="vertical",
                        aspect=40, extend='both',)
    cbar.set_label('% Over Background', fontsize='x-large')

    return fig


def fig2(ds,
         vlims,
         lons,
         at_alt=None,
         tlims=None,
         numplots=6,
         gitm=True,
         mag_eq=None,
         ):

    fig, axs = plt.subplots(3, 2, figsize=(15, 10),
                            sharex=True,
                            sharey=True,)

    # plotting...
    # Setup some things...
    t = 0

    axes = axs.flatten()

    if tlims is not None:
        ds = ds.sel(time=slice(tlims[0], tlims[1]))

    while t < numplots:
        added_lon = False
        if ds.sel(lon=lons[t], method='nearest').sum() == 0:
            lons[t] += 4
            added_lon = True
        im = ds.sel(lon=lons[t], method='nearest').plot(x='time',
                                                        vmin=-vlims,
                                                        vmax=vlims,
                                                        ax=axes[t],
                                                        add_colorbar=False,
                                                        cmap='bwr'
                                                        )
        if added_lon:
            lons[t] -=4

        if lons[t] < 180 and lons[t] > 0:
            axes[t].set_title("Geographic Longitude = %i° East" %
                          int(lons[t]), fontsize='x-large')
        elif lons[t] > 180:
            axes[t].set_title("Geographic Longitude = %i° West" %
                          (360-int(lons[t])), fontsize='x-large')
        else:
            axes[t].set_title("Geographic Longitude = %i°" %
                          int(lons[t]), fontsize='x-large')
            
        axes[t].xaxis.set_major_formatter(DateFormatter('%H:%M'))
        
        axes[t].vlines(pd.Timestamp('2011-05-21 12:00:00'),-90,90,
                       color='peru', linestyle = 'solid', linewidth=3)
        
        lts = ut_to_lt(ds.time.values, lons[t])
        dawn = np.argmin(np.abs(lts-6))
        dusk = np.argmin(np.abs(lts-18))
        if dawn > 2 and dawn < len(ds.time.values)-2:
            axes[t].vlines(ds.time.values[dawn],-90,90,
                      color='lime', linestyle = 'dotted', linewidth=3)
        if dusk > 2 and dusk < len(ds.time.values)-2:
            axes[t].vlines(ds.time.values[dusk],-90,90,
                      color='lime', linestyle = 'dashed', linewidth=3)
            
        if mag_eq is not None:
            axes[t].hlines(
                mag_eq['mlat'][np.argmin(np.abs(mag_eq['glon'] - lons[t]))],
                ds.time.values[0], ds.time.values[-1],
                            linestyle='--', color='k')
            
        # axes[t].axvspan(ds.time.values[dawn], ds.time.values[dawn], alpha=0.15, color='k')

        axes[t].set_ylabel('')
        axes[t].set_xlabel('')
        t += 1

    fig.supxlabel('Universal Time (HH:MM)', fontsize='x-large')
    fig.supylabel('Latitude (Degrees North)', fontsize='x-large')

    if gitm:
        fig.suptitle("Perturbation Total Neutral Mass Density from GITM at %ikm" % at_alt,
                 fontsize='xx-large', fontweight='heavy')
    else:
        fig.suptitle("Perturbation TEC from SAMI",
                 fontsize='xx-large', fontweight='heavy')

    fig.tight_layout()

    cbar = fig.colorbar(im, ax=axs.ravel().tolist(), orientation="vertical",
                        aspect=40, extend='both',)
    if gitm:
        cbar.set_label('% Over Background', fontsize='x-large')
    else:
        cbar.set_label('Perturbation TEC (TECUnits)', fontsize='x-large')
    return fig


def fig3(da,
         vmin,
         vmax,
         suptitle,
         dial_cmap='jet',
         dial_kwargs={},
         ):

    # for local time:
    time_here = pd.Timestamp(da.time.values)
    # Find central longitude (midnight Local time)
    lons = np.arange(0, 360)
    lts = ut_to_lt([time_here], lons)
    central_lon = lons[np.argmin(np.abs(24-lts))]

    fig = plt.figure(figsize=(10, 4))
    ax1 = fig.add_subplot(1, 2, 1,
                          projection=ccrs.Orthographic(central_lon, 90))
    da.plot(ax=ax1, x='lon', transform=ccrs.PlateCarree(),
                   cmap=dial_cmap,
                   vmin=vmin, vmax=vmax, add_colorbar=False,
            **dial_kwargs)
    ax1.set_title('Northern Hemisphere')
    ax1.add_feature(Nightshade(time_here), alpha=0.3)
    ax1.gridlines(color='black', linestyle='--',
                              alpha=0.6, draw_labels=True)
    ax1.coastlines(zorder=3, color='black', alpha=1)


    ax2 = fig.add_subplot(1, 2, 2,
                          projection=ccrs.Orthographic(
                              central_lon-180, -90))
    im = da.plot(ax=ax2, x='lon', transform=ccrs.PlateCarree(),
                   cmap=dial_cmap,
            add_colorbar=False, **dial_kwargs)
    ax2.set_title('Southern Hemisphere')
    ax2.add_feature(Nightshade(time_here), alpha=0.3)
    ax2.gridlines(color='black', linestyle='--',
                              alpha=0.6, draw_labels=True)
    ax2.coastlines(zorder=3, color='black', alpha=1)
    fig.suptitle(suptitle, fontsize='xx-large', fontweight='heavy')

    fig.tight_layout()

    cbar = fig.colorbar(im, ax=[ax1, ax2], orientation="vertical",
                        aspect=17, extend='max')
    cbar.set_label(
        'Altitude Integrated \nJoule Heating ($\mathrm{W/M^2}$)', fontsize='x-large')

    return fig


def fig4(ds,
         vlims,
         lons,
         tlims=None,
         numplots=6,
         gitm=False,
         mag_eq=None
         ):
    
    fig = fig2(ds,
         vlims,
         lons,
         tlims=tlims,
         numplots=numplots,
         gitm=False,
              mag_eq=mag_eq)
    return fig


def fig5(ds,
         vlims,
         times,
         storm_start,
         lon,
         numplots=6,
         ):

    fig, axs = plt.subplots(3, 2, figsize=(13, 10),
                            sharex=True,
                            sharey=True,)

    # plotting...
    # Setup some things...
    t = 0

    axes = axs.flatten()

    while t < numplots:
        im = ds.sel(time=times[t], method='nearest').interp(
            lat=np.arange(-90, 90, 1),
            alt=np.arange(150, 1500, 20)).plot(x='lat',
                                                        vmin=-vlims,
                                                        vmax=vlims,
                                                        ax=axes[t],
                                                        add_colorbar=False,
                                                        cmap='bwr'
                                                        )

        axes[t].set_title("%s From Storm Onset (%s Local Time)"
                          %(timedeltatotime(times[t]-storm_start),
                            ut_to_lt([times[t]], lon)[0]), fontsize = 'x-large')

        axes[t].set_ylabel('')
        axes[t].set_xlabel('')
        axes[t].set_ylim(100, 1000)
        t += 1

    fig.supxlabel('Latitude (Degrees North)', fontsize = 'x-large')
    fig.supylabel('Altitude (km)', fontsize = 'x-large')

    fig.suptitle("Perturbation Electron Density from SAMI\n"
                 "Along %s Geographic Longitude"
                 %(str(int(lon)) + "° East" if lon < 180 else str(360-int(lon)) + "° West"),
                 fontsize = 'xx-large', fontweight = 'heavy')

    fig.tight_layout()

    cbar=fig.colorbar(im, ax = axs.ravel().tolist(), orientation = "vertical",
                        aspect = 40, extend = 'both',)
    cbar.set_label('% Over Background', fontsize = 'x-large')

    return fig