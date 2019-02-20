#include <fstream>
#include <iostream>
#include <iomanip>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_LevelBld.H>

using namespace amrex;

#include <Transport_F.H>
#include <main_F.H>
#include <PlotFileFromMF.H>
#include <actual_Creactor.h>

/**********************************/
int main (int argc, 
          char* argv[])
{
    amrex::Initialize(argc,argv);

    Real strt_time, stop_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();

    std::cout << std::setprecision(15);

    int n_cell, max_grid_size;
    int cvode_meth,cvode_itmeth,cvode_iJac,cvode_iE,cvode_iDense,ivode;
    int write_plotfile;
    Real time, dt;
    int ndt;
    bool do_tiling;
    std::string probin_file, pltfile;
    std::string txtfile_in=""; 

    // inputs parameters
    {
      // ParmParse is way of reading inputs from the inputs file
      ParmParse pp;
      
      // probin file
      pp.get("probin_file",probin_file);

      // We need to get n_cell from the inputs file - this is the number of
      // cells on each side of a square (or cubic) domain.
      pp.get("n_cell",n_cell);

      // Default nsteps to 0, allow us to set it to something else in the
      // inputs file
      pp.get("max_grid_size",max_grid_size);

      // Total nb of integration dt to take
      pp.get("ndt",ndt);
      // Time steps
      pp.get("dt",dt);
      // Select CVODE solve method.
      //   1 = Adams (for non-stiff problems)
      //   2 = BDF (for stiff problems)
      pp.get("cvode_meth",cvode_meth);
      // Select CVODE solver iteration method.
      //   1 = Functional iteration
      //   2 = Newton iteration
      pp.get("cvode_itmeth",cvode_itmeth);

      cvode_iJac = 0;
      pp.query("cvode_iJac",cvode_iJac);
      // Select CVODE Jacobian eval.
      //   0 = finite differences
      //   1 = User supplied function
      
      pp.get("cvode_iE",cvode_iE);
      // Select CVODE type of energy employed.
      //1 for UV, 2 for HP
      //   1 = Internal energy
      //   anything else = enthalpy (PeleLM restart)
      
      pp.get("cvode_iDense",cvode_iDense);
      //1 for regular dense solver, otherwise it is sparse
      
      pp.get("ivode",ivode);
      // Select bet CVODE (=1) or dvode (=other)

      pp.get("write_plotfile",write_plotfile);
      pp.get("do_tiling",do_tiling);

      ParmParse ppa("amr");
      ppa.query("plot_file",pltfile); 

      pp.query("txtfile_in",txtfile_in);

    }

    // Print parameters on screen
    if (cvode_meth < 1)
      amrex::Abort("Unknown cvode_meth");
    if (cvode_itmeth < 1)
      amrex::Abort("Unknown cvode_itmeth");

    amrex::Print() << "Problem domain size: nx = ny = nz = " << n_cell << std::endl;
    amrex::Print() << "Max grid size: " << max_grid_size << std::endl;
    //amrex::Print() << "Total integration time: " << dt << std::endl;
    amrex::Print() << "External time-stepping: " << dt << std::endl;
    amrex::Print() << "Final time: " << ndt*dt << std::endl;
    amrex::Print() << "CVODE method (if dvode is used stiff method employed): ";
    if (cvode_meth == 1) {
      amrex::Print() << "Adams (non-stiff)";
    } else if (cvode_meth == 2) {
        amrex::Print() << "BDF (stiff)";
    }
    amrex::Print() << std::endl;
    amrex::Print() << "CVODE iteration method (if dvode is used Newton is employed): ";
    if (cvode_itmeth == 1) {
      amrex::Print() << "Functional";
    } else if (cvode_itmeth == 2) {
        amrex::Print() << "Newton";
    }
    amrex::Print() << std::endl;
    amrex::Print() << "CVODE use of Analytical J (disabled if dvode employed): ";
    if (cvode_iJac == 0) {
      amrex::Print() << "NO";
    } else {
        amrex::Print() << "YUP";
    }
    
    amrex::Print() << std::endl;


    // take care of probin init to initialize problem
    int probin_file_length = probin_file.length();
    std::vector<int> probin_file_name(probin_file_length);
    for (int i = 0; i < probin_file_length; i++)
	    probin_file_name[i] = probin_file[i];
    extern_init(&(probin_file_name[0]),&probin_file_length, &cvode_iE);
    extern_cInit(&cvode_meth, &cvode_itmeth,&cvode_iJac, &cvode_iE, &cvode_iDense, &max_grid_size);


    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    IntVect dom_lo(IntVect(D_DECL(0,0,0)));
    //IntVect dom_hi(IntVect(D_DECL(n_cell-1, 0, 0)));
    IntVect dom_hi(IntVect(D_DECL(n_cell-1, 4, 0)));
    Box domain(dom_lo, dom_hi);

    // Initialize the boxarray "ba" from the single box "bx"
    ba.define(domain);

    // Break up boxarray "ba" into chunks no larger than "max_grid_size"
    // along a direction
    const IntVect ChunkSize = IntVect::TheUnitVector(); //parent->maxGridSize(level);
    IntVect chunk(ChunkSize);
    chunk[0] = 1; //max_grid_size;
    chunk[1] = max_grid_size; 
    chunk[2] = 1; 
    ba.maxSize(chunk);

    // This defines the physical size of the box.  Right now the box is
    // [-1,1] in each direction.
    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
      real_box.setLo(n,-1.0);
      real_box.setHi(n, 1.0);
    }

    // This sets the boundary conditions to be doubly or triply periodic
    int is_periodic[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) {
      is_periodic[i] = 1;
    }

    // This defines a Geometry object
    geom.define(domain,&real_box,CoordSys::cartesian,is_periodic);

    // This is a bit useless: dx is cst
    //std::vector<Real> plo(3,0), phi(3,0), dx(3,1);
    Real plo ;
    //for (int i=0; i<BL_SPACEDIM; ++i) {
    //        phi[i] = domain.length(i); 
    //        dx[i] = (phi[i] - plo[i])/domain.length(i); 
    //  }

    // Ncomp = number of species to integrate
    int Ncomp;
    get_num_spec(&Ncomp);

    // time = starting time in the simulation
    time = 0.0;

    // Create MultiFabs with no ghost cells.
    DistributionMapping dm(ba);
    MultiFab mf(ba, dm, Ncomp+1, 0);
    MultiFab rY_source_ext(ba,dm,Ncomp,0);
    MultiFab mfE(ba, dm, 1, 0);
    MultiFab rY_source_energy_ext(ba,dm,1,0);
    MultiFab temperature(ba,dm,1,0);
    iMultiFab mask(ba,dm,1,0);
    MultiFab cost(ba,dm,1,0);

    mask.setVal(1,0,1,0);

    amrex::Print() << " " << std::endl;
    amrex::Print() << "... Top chrono ..." << std::endl;
    amrex::Print() << " " << std::endl;
    strt_time = ParallelDescriptor::second();

    // This is first option: initialize via txt file 
    // obtained either via chemistry solver or fextract
    // !! Check data format in main_nd.F read_data_from_txt
    if (txtfile_in != "")
    {
         int txtfile_in_length = txtfile_in.length();
         std::vector<int> txtfile_in_name(txtfile_in_length);
         for (int i = 0; i < txtfile_in_length; i++)
	 txtfile_in_name[i] = txtfile_in[i];

	 printf(" (main) Will read data_from_txt ... \n");
	 read_data_from_txt(&(txtfile_in_name[0]),&txtfile_in_length,&plo);

#ifdef _OPENMP
#pragma omp parallel
#endif
	 // Initialize MF with data from txt file
	 // Note rY_source_ext is now in Fortran module so could technically
	 // be removed from here
         int count = 0; 
         for ( MFIter mfi(mf, do_tiling); mfi.isValid(); ++mfi )
            {
              const Box& box = mfi.tilebox();

              initialize_data(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
                		BL_TO_FORTRAN_N_3D(mf[mfi],0),
                		BL_TO_FORTRAN_N_3D(rY_source_ext[mfi],0),
			        BL_TO_FORTRAN_N_3D(mfE[mfi],0),
			        BL_TO_FORTRAN_N_3D(rY_source_energy_ext[mfi],0));
              count = count + 1;
            }

         if (count != 1) 
         {
              amrex::Print() << "INFO: "<< count <<" FabArray in MF ... "<< std::endl;
         }
    } else {

#ifdef _OPENMP
#pragma omp parallel
#endif
        //int count = 1; 
        for ( MFIter mfi(mf, do_tiling); mfi.isValid(); ++mfi )
        {
          //amrex::Print() << "... working on mf number " << count << std::endl;
          const Box& box = mfi.tilebox();

          initialize_data_byhand(ARLIM_3D(box.loVect()), ARLIM_3D(box.hiVect()),
            		BL_TO_FORTRAN_N_3D(mf[mfi],0),
			BL_TO_FORTRAN_N_3D(rY_source_ext[mfi],0), 
			BL_TO_FORTRAN_N_3D(mfE[mfi],0),
			BL_TO_FORTRAN_N_3D(rY_source_energy_ext[mfi],0),
			&plo);
          //count = count + 1;
        }

    }


    // Check initialisation
    if (write_plotfile)
    {
      std::string outfile = Concatenate(pltfile,1); // Need a number other than zero for reg test to pass
      PlotFileFromMF(mf,outfile);
    }


#ifdef _OPENMP
#pragma omp parallel
#endif
    /* Advance the MF by dt*ndt */
    std::ofstream myfile;
    myfile.open ("fort.txt");
    myfile.precision(12);
    int count = 1;
    int count_box;
    for ( MFIter mfi(mf, do_tiling); mfi.isValid(); ++mfi )
    {
        amrex::Print() << "\n ... working on mf number " << count << std::endl;

        const Box& box = mfi.tilebox();
	FArrayBox& Fb = mf[mfi];
	FArrayBox& Fbsrc = rY_source_ext[mfi];
	FArrayBox& FbE = mfE[mfi];
	FArrayBox& FbEsrc = rY_source_energy_ext[mfi];

        /* Pack the data */
	count_box = 1;
	// rhoY,T
	double tmp_vect[max_grid_size*(Ncomp+1)];
	// rhoY_src_ext
	double tmp_src_vect[max_grid_size*(Ncomp)];
	// rhoE/rhoH
	double tmp_vect_energy[max_grid_size];
	double tmp_src_vect_energy[max_grid_size];
	for (BoxIterator bit(box); bit.ok(); ++bit) {
		tmp_vect_energy[(count_box-1)] = FbE(bit(),0);
		tmp_src_vect_energy[(count_box-1)] = FbEsrc(bit(),0);
		for (int i=0;i<Ncomp; i++){
			tmp_vect[(count_box-1)*(Ncomp+1) + i] = Fb(bit(),i);
			tmp_src_vect[(count_box-1)*(Ncomp) + i] = Fbsrc(bit(),i);
		}
	        tmp_vect[(count_box-1)*(Ncomp+1) + Ncomp] = Fb(bit(), Ncomp);
		count_box = count_box + 1;
	}
        amrex::Print() << "... will contain " << count_box-1 << "cells" << std::endl;

        /* Solve the problem */
	Real time_tmp, dt_incr;
	dt_incr =  dt;
	time_tmp = time;
	int reInit = 1;
	//printf("#TIME TEMPERATURE \n");
	myfile << "#TIME TEMPERATURE P Yks \n";
	for (int i = 0; i < ndt; ++i) {
		actual_cReact(tmp_vect, tmp_src_vect, 
				tmp_vect_energy, tmp_src_vect_energy,
				&plo, &dt_incr, &time_tmp, &reInit);
				//&plo, &dt_incr, &time, &reInit);
	        // increment time with true dt_incr
		time_tmp = time_tmp + dt_incr;
		// fix new dt_incr to chosen value, hoping cvode will reach it
		dt_incr = dt;
                //printf("--> at time : %2.8e, ",time_tmp);
		//printf(" out state (T) is %4.4f and Ysrc(OH) is \n", tmp_vect[Ncomp], plo);
		//printf("%2.8e %4.4f \n", time_tmp, tmp_vect[Ncomp]);
		myfile <<  time_tmp << " " << tmp_vect[Ncomp] << " " << plo <<" ";
		for (int ij = 0; ij < Ncomp; ++ij) {
			myfile << tmp_vect[ij]<< " ";
		}
		myfile <<"\n";
	}
        count = count+1;

        /* Unpack the data ? */
	count_box = 1;
	for (BoxIterator bit(box); bit.ok(); ++bit) {
		for (int i=0;i<Ncomp+1; i++){
			Fb(bit(),i) = tmp_vect[(count_box-1)*(Ncomp+1) + i];
		}
		count_box = count_box + 1;
	}
    }
    myfile.close();

    stop_time = ParallelDescriptor::second();
    ParallelDescriptor::ReduceRealMax(stop_time,IOProc);
    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Time = " << stop_time - strt_time << std::endl;

    //MultiFab::Copy(temperature,mf,Ncomp+1,0,1,0);
    // Check results
    if (write_plotfile)
    {
      std::string outfile = Concatenate(pltfile,2); // Need a number other than zero for reg test to pass
      PlotFileFromMF(mf,outfile);
      //PlotFileFromMF(temperature,outfile);
    }

    //printf("THER");
    extern_close();
    extern_cFree();
    amrex::Finalize();
    return(0);
}
