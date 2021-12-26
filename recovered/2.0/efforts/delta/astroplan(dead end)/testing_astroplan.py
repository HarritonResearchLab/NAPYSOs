def first_plot():
    #Import(s)
    from astroplan.plots import plot_finder_image
    from astroplan import FixedTarget
    import matplotlib.pyplot as plt

    #Action

    ngc7000 = FixedTarget.from_name('NGC 7000')
    ax, hdu = plot_finder_image(ngc7000)
    plt.show()

first_plot()