package main

import (
	"fmt"
	"math"
	"time"

	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
)

func main() {
	nu := 0.5
	dt := 0.01
	T := 500
	_ = T
	L := 2 * math.Pi
	num_grids := math.Ceil(L / (dt / nu))
	dx := L / float64(num_grids)

	var x [315]float64
	var val float64 = 0

	for i := 0; i < int(num_grids); i++ {
		x[i] = val
		val += dx
	}

	fmt.Println("Step size =", dx, "s")
	nx := len(x)
	fmt.Println("Number of x nodes =", nx)

	var u0 [315]float64

	for i := 0; i < int(num_grids); i++ {
		u0[i] = math.Sin(x[i]) + 0.5*math.Sin(0.5*x[i])
	}

	u := u0
	var ustar [315]float64
	var unp1 [315]float64

	var mac_cor_u_0, mac_cor_u_1, mac_cor_u_3, mac_cor_u_5, mac_cor_time = mac_cormack(T, nx, ustar, u, unp1, nu)
	var lax_u_0, lax_u_1, lax_u_3, lax_u_5, lax_time = lax(T, nx, ustar, u, unp1, nu)
	var lax_wend_u_0, lax_wend_u_1, lax_wend_u_3, lax_wend_u_5, lax_wend_time = lax_wend(T, nx, ustar, u, unp1, nu)

	graph_plt("MacCormack Method", x, mac_cor_u_0, mac_cor_u_1, mac_cor_u_3, mac_cor_u_5)
	graph_plt("Lax Method", x, lax_u_0, lax_u_1, lax_u_3, lax_u_5)
	graph_plt("Lax Wendroff Method", x, lax_wend_u_0, lax_wend_u_1, lax_wend_u_3, lax_wend_u_5)

	plt2 := plot.New()

	plt2.Title.Text = "Solution at t = 5s"
	plt2.X.Label.Text = "X"
	plt2.Y.Label.Text = "U"
	plt2.Y.Max = 1.5
	plt2.X.Max = 7

	err := plotutil.AddLinePoints(plt2,
		"u", put_in_plots(x[:], u[:]),
		"MacCormack", put_in_plots(x[:], mac_cor_u_5[:]),
		"Lax-Wendroff", put_in_plots(x[:], lax_wend_u_5[:]),
		"Lax", put_in_plots(x[:], lax_u_5[:]))

	if err != nil {
		panic(err)
	}
	// Save the plot to a PNG file.
	if err := plt2.Save(10*vg.Inch, 6*vg.Inch, "solat5.png"); err != nil {
		panic(err)
	}
	fmt.Println("Mac Cormack Time:", mac_cor_time)
	fmt.Println("Lax Method Time:", lax_time)
	fmt.Println("Lax Wendroff Time:", lax_wend_time)
}

func put_in_plots(x_input []float64, u_input []float64) plotter.XYs {
	pts := make(plotter.XYs, 315)
	for i := range pts {
		pts[i].X = x_input[i]
		pts[i].Y = u_input[i]
	}
	return pts
}

func graph_plt(graph_title string, x [315]float64, u_0 [315]float64, u_1 [315]float64, u_3 [315]float64, u_5 [315]float64) {
	plt := plot.New()

	plt.Title.Text = graph_title
	plt.X.Label.Text = "X"
	plt.Y.Label.Text = "U"
	plt.Y.Max = 1.5
	plt.X.Max = 7

	err := plotutil.AddLinePoints(plt,
		"t = 0 s", put_in_plots(x[:], u_0[:]),
		"t = 1 s", put_in_plots(x[:], u_1[:]),
		"t = 3 s", put_in_plots(x[:], u_3[:]),
		"t = 5 s", put_in_plots(x[:], u_5[:]))

	if err != nil {
		panic(err)
	}
	graph_title = graph_title + ".png"
	// Save the plot to a PNG file.
	if err := plt.Save(10*vg.Inch, 8*vg.Inch, graph_title); err != nil {
		panic(err)
	}
}

func mac_cormack(T int, nx int, ustar [315]float64, u [315]float64, unp1 [315]float64, nu float64) (u_0 [315]float64, u_1 [315]float64, u_3 [315]float64, u_5 [315]float64, elapsed time.Duration) {
	start := time.Now()
	for t := 0; t < T; t++ {

		if t%2 == 0 {
			for i := 0; i < nx-1; i++ {
				ustar[i] = u[i] - 0.5*nu*(math.Pow(u[i+1], 2)-math.Pow(u[i], 2))
			}
			ustar[313] = u[313] - 0.5*nu*(math.Pow(u[1], 2)-math.Pow(u[313], 2))

			unp1[0] = 0.5 * (u[0] + ustar[0] - 0.5*nu*(math.Pow(ustar[0], 2)-math.Pow(ustar[nx-2], 2)))

			for i := 1; i < nx; i++ {
				unp1[i] = 0.5 * (u[i] + ustar[i] - 0.5*nu*(math.Pow(ustar[i], 2)-math.Pow(ustar[i-1], 2)))
			}
		} else {

			ustar[0] = u[0] - 0.5*nu*(math.Pow(u[0], 2)-math.Pow(u[nx-2], 2))

			for i := 1; i < nx; i++ {
				ustar[i] = u[i] - 0.5*nu*(math.Pow(u[i], 2)-math.Pow(u[i-1], 2))
			}

			for i := 0; i < nx-1; i++ {
				unp1[i] = 0.5 * (u[i] + ustar[i] - 0.5*nu*(math.Pow(ustar[i+1], 2)-math.Pow(ustar[i], 2)))
			}

			unp1[313] = 0.5 * (u[313] + ustar[313] - 0.5*nu*(math.Pow(ustar[1], 2)-math.Pow(ustar[313], 2)))
		}

		u = unp1

		if t%100 == 0 {
			if float32(t)*0.01 == 0. {
				u_0 = u
			} else if float32(t)*0.01 == 1. {
				u_1 = u
			} else if float32(t)*0.01 == 3. {
				u_3 = u
			}

		}

	}
	u_5 = u
	elapsed = time.Since(start)
	return

}

func lax_wend(T int, nx int, ustar [315]float64, u [315]float64, unp1 [315]float64, nu float64) (u_0 [315]float64, u_1 [315]float64, u_3 [315]float64, u_5 [315]float64, elapsed time.Duration) {
	start := time.Now()
	for t := 0; t < T; t++ {
		Ajp := 0.5 * (u[0] + u[1])
		Ajm := 0.5 * (u[0] + u[nx-2])
		unp1[0] = u[0] - 0.25*nu*(math.Pow(u[1], 2)-math.Pow(u[nx-2], 2)) + 0.5*math.Pow(nu, 2)*(0.5*Ajp*(math.Pow(u[1], 2)-math.Pow(u[0], 2))-0.5*Ajm*(math.Pow(u[0], 2)-math.Pow(u[nx-2], 2)))
		for i := 1; i < nx-1; i++ {
			Ajp = 0.5 * (u[i] + u[i+1])
			Ajm = 0.5 * (u[i] + u[i-1])
			unp1[i] = u[i] - 0.25*nu*(math.Pow(u[i+1], 2)-math.Pow(u[i-1], 2)) + 0.5*math.Pow(nu, 2)*(0.5*Ajp*(math.Pow(u[i+1], 2)-math.Pow(u[i], 2))-0.5*Ajm*(math.Pow(u[i], 2)-math.Pow(u[i-1], 2)))
		}

		Ajp = 0.5 * (u[313] + u[1])
		Ajm = 0.5 * (u[313] + u[313-1])
		unp1[313] = u[313] - 0.25*nu*(math.Pow(u[1], 2)-math.Pow(u[313-1], 2)) + 0.5*math.Pow(nu, 2)*(0.5*Ajp*(math.Pow(u[1], 2)-math.Pow(u[313], 2))-0.5*Ajm*(math.Pow(u[313], 2)-math.Pow(u[313-1], 2)))
		u = unp1

		if t%100 == 0 {
			if float32(t)*0.01 == 0. {
				u_0 = u
			} else if float32(t)*0.01 == 1. {
				u_1 = u
			} else if float32(t)*0.01 == 3. {
				u_3 = u
			}

		}
	}

	u_5 = u
	elapsed = time.Since(start)
	return
}

func lax(T int, nx int, ustar [315]float64, u [315]float64, unp1 [315]float64, nu float64) (u_0 [315]float64, u_1 [315]float64, u_3 [315]float64, u_5 [315]float64, elapsed time.Duration) {
	start := time.Now()
	for t := 0; t < T; t++ {
		unp1[0] = 0.5*(u[1]+u[nx-2]) - 0.25*nu*(math.Pow(u[1], 2)-math.Pow(u[nx-2], 2))
		for i := 1; i < nx-1; i++ {
			unp1[i] = 0.5*(u[i+1]+u[i-1]) - 0.25*nu*(math.Pow(u[i+1], 2)-math.Pow(u[i-1], 2))
		}
		unp1[313] = 0.5*(u[1]+u[313-1]) - 0.25*nu*(math.Pow(u[1], 2)-math.Pow(u[313-1], 2))
		u = unp1

		if t%100 == 0 {
			if float32(t)*0.01 == 0. {
				u_0 = u
			} else if float32(t)*0.01 == 1. {
				u_1 = u
			} else if float32(t)*0.01 == 3. {
				u_3 = u
			}

		}

	}
	u_5 = u
	elapsed = time.Since(start)
	return
}
