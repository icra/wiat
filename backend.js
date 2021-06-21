/*
  Draft for WIAT backend
*/

/*
  UTILS
*/
  //sum array of numbers
  Array.prototype.sum=function(){return this.reduce((p,c)=>(p+c),0)};

/*
  ASSESSMENT: main class
*/
class Assessment{
  constructor(){
    this.name="new assessment";

    //assessment period
    this.assessment_period_start = new Date();
    this.assessment_period_end   = new Date();

    //array for "Industry" objects
    this.industries=[];
  }

  /*
    methods
  */

  //total emissions (all industries)
  TotalGHG(){
    return this.industries.map(i=>i.wwt_KPI_GHG().total).sum(); //kgCO2eq
  }

  //total energy consumed (all industries)
  TotalNRG(){
    return this.industries.map(i=>i.wwt_nrg_cons).sum(); //kWh
  }

  //add industry
  add_industry(){
    let i = new Industry();
    this.industries.push(i);
  }

  //delete industry
  delete_industry(i){
    this.industries.splice(i,1);
  }
};

/*
  INDUSTRY: industries live inside an assessment object
*/
class Industry{
  constructor(){
    this.name="new industry";

    //energy
    this.wwt_nrg_cons = 0; //kWh | energy consumed from the grid
    this.wwt_conv_kwh = 0; //kgCO2eq/kWh | conversion factor

    //BOD (creates CH4)
    this.wwt_bod_infl = 0; //kgBOD
    this.wwt_bod_effl = 0; //kgBOD

    //TN (creates N2O)
    this.wwt_tn_infl = 0; //kgN
    this.wwt_tn_effl = 0; //kgN

    //emission factors (treatment)
    this.wwt_ch4_efac_tre = 0; //kgCH4/kgBOD
    this.wwt_n2o_efac_tre = 0; //kgN2O-N/kgN

    //emission factors (discharge)
    this.wwt_ch4_efac_dis = 0; //kgCH4/kgBOD
    this.wwt_n2o_efac_dis = 0; //kgN2O-N/kgN

    //fuel engines
    this.wwt_vol_fuel = 0; //L of fuel
    this.wwt_fuel_typ = 0; //Option | type of fuel (see Tables)

    //biogas
    this.wwt_biog_pro = 0;  //Nm3 | total biogas produced
    this.wwt_biog_fla = 98; //% of biogas produced that is flared
    this.wwt_biog_val = 0;  //% of biogas produced that is used for heat
    this.wwt_biog_lkd = 2;  //% of biogas produced that is leaked
    this.wwt_ch4_biog = 59; //% of CH4 in biogas (volume)
    this.wwt_dige_typ = 0;  //Option | type of fuel for digester
    this.wwt_fuel_dig = 0;  //L | volume of fuel used in the digester

    //fuel used in water reuse trucks
    this.wwt_reus_trck_typ = 0; //Option | type of fuel
    this.wwt_reus_vol_trck = 0; //L | volume of fuel used

    //SLUDGE MANAGEMENT
    this.wwt_mass_slu = 0; //kg | raw sludge removed from wwtp as dry mass
    this.wwt_bod_slud = 0; //kg | BOD removed as sludge

    //sludge storage
    this.wwt_mass_slu_sto  = 0; //kg of sludge stored
    this.wwt_time_slu_sto  = 0; //days
    this.wwt_slu_sto_TVS   = 0; //%
    this.wwt_slu_sto_f_CH4 = 0; //% for CH4 potential
    this.wwt_slu_sto_EF    = 0; //%

    //sludge composting
    this.wwt_mass_slu_comp                          = 0; //kg of sludge composted
    this.wwt_slu_comp_emis_treated_or_piles_covered = 0; //yes/no
    this.wwt_slu_comp_solids_content                = 0; //%
    this.wwt_slu_comp_TVS                           = 0; //%
    this.wwt_slu_comp_N_cont                        = 0; //%
    this.wwt_slu_comp_low_CN_EF                     = 0.015; //kgN2O-N/kgN
    this.wwt_slu_comp_uncovered_pile_EF             = 0.025; //kgCH4/kgC
    this.wwt_slu_comp_seqst_rate                    = 0.25; //kgCO2eq/kgSludge

    //sludge incineration
    this.wwt_mass_slu_inc   = 0;    //kg of sludge incinerated
    this.wwt_temp_inc       = 1023; //K | temperature incineration
    this.wwt_slu_inc_N_cont = 0;    //% of N
    this.wwt_slu_inc_SNCR   = 0;    //boolean

    //sludge LA
    this.wwt_mass_slu_app          = 0; //kg of sludge sent to LA
    this.wwt_slu_la_solids_content = 0; //%
    this.wwt_slu_la_TVS            = 0; //%
    this.wwt_slu_la_N_cont         = 0; //%
    this.wwt_slu_la_EF             = 0; //gN2O-N/gN

    //sludge LF
    this.wwt_mass_slu_land      = 0;    //kg of sludge sent to LF
    this.wwt_slu_lf_TVS         = 0;    //%
    this.wwt_slu_lf_uncertainty = 0.9;  //adimensional
    this.wwt_slu_lf_CH4_in_gas  = 50;   //%
    this.wwt_slu_lf_DOCf        = 80;   //%
    this.wwt_slu_lf_decomp_3yr  = 69.9; //%
    this.wwt_slu_lf_MCF         = 1;    //ratio
    this.wwt_slu_lf_N_cont      = 0;    //% N content
    this.wwt_slu_lf_low_CN_EF   = 0.015; //kgN2O-N/kgN

    //sludge SP
    this.wwt_mass_slu_stock = 0;  //kg of sludge stockpiled
    this.wwt_slu_sp_lifespan = 0; //years

    //sludge truck transport
    this.wwt_trck_typ = 0; //Option | fuel type
    this.wwt_vol_tslu = 0; //L | volume of fuel
  }

  /*
    GHG emissions (kgCO2eq)
  */
    //total GHG emissions
    wwt_KPI_GHG(){
      let sources=[
        this.wwt_KPI_GHG_elec(),
        this.wwt_KPI_GHG_fuel(),
        this.wwt_KPI_GHG_tre(),
        this.wwt_KPI_GHG_biog(),
        this.wwt_KPI_GHG_slu(),
        this.wwt_KPI_GHG_reus_trck(),
        this.wwt_KPI_GHG_disc(),
      ];

      //gases (numbers)
      let co2 = sources.map(s=>s.co2).sum();
      let ch4 = sources.map(s=>s.ch4).sum();
      let n2o = sources.map(s=>s.n2o).sum();

      //total
      let total = sources.map(s=>s.total).sum();
      return {total,co2,ch4,n2o};
    }

    //indirect emissions from electricity consumption
    wwt_KPI_GHG_elec(){
      let co2 = this.wwt_nrg_cons*this.wwt_conv_kwh;
      let ch4 = 0;
      let n2o = 0;
      let total = co2+ch4+n2o;
      return {total,co2,ch4,n2o};
    }

    //emissions from fuel engines
    wwt_KPI_GHG_fuel(){
      let vol   = this.wwt_vol_fuel;
      let fuel  = Tables.get_row('Fuel type',this.wwt_fuel_typ); //object
      let co2   = vol*fuel.FD*fuel.NCV/1000*fuel.EFCO2;
      let ch4   = vol*fuel.FD*fuel.NCV/1000*fuel.EFCH4.engines*Cts.ct_ch4_eq.value;
      let n2o   = vol*fuel.FD*fuel.NCV/1000*fuel.EFN2O.engines*Cts.ct_n2o_eq.value;
      let total = co2+n2o+ch4;
      return {total,co2,ch4,n2o};
    }

    //emissions from biogas (fuel used in digester)
    wwt_KPI_GHG_biog_dig(){
      let vol   = this.wwt_fuel_dig;
      let fuel  = Tables.get_row('Fuel type',this.wwt_dige_typ);
      let co2   = vol*fuel.FD*fuel.NCV/1000*fuel.EFCO2
      let ch4   = vol*fuel.FD*fuel.NCV/1000*fuel.EFCH4.engines*Cts.ct_ch4_eq.value;
      let n2o   = vol*fuel.FD*fuel.NCV/1000*fuel.EFN2O.engines*Cts.ct_n2o_eq.value;
      let total = co2+n2o+ch4;
      return {total,co2,ch4,n2o};
    }

    //emissions from treatment
    wwt_KPI_GHG_tre(){
      let co2   = 0;
      let ch4   = (this.wwt_bod_infl-this.wwt_bod_slud)*this.wwt_ch4_efac_tre*Cts.ct_ch4_eq.value;
      let n2o   = this.wwt_tn_infl*this.wwt_n2o_efac_tre*Cts.ct_N_to_N2O_44_28.value*Cts.ct_n2o_eq.value;
      let total = co2+ch4+n2o;
      return {total,co2,ch4,n2o};
    }

    //emissions from biogas
    wwt_KPI_GHG_biog(){
      let sources=[
        this.wwt_KPI_GHG_biog_flared(),
        this.wwt_KPI_GHG_biog_valorized(),
        this.wwt_KPI_GHG_biog_leaked(),
        this.wwt_KPI_GHG_biog_dig(),
      ];

      //gases (numbers)
      let co2 = sources.map(s=>s.co2).sum();
      let ch4 = sources.map(s=>s.ch4).sum();
      let n2o = sources.map(s=>s.n2o).sum();

      //total
      let total = sources.map(s=>s.total).sum();
      return {total,co2,ch4,n2o};
    }

    //emissions from biogas flared
    wwt_KPI_GHG_biog_flared(){
      let moles_biogas        = this.wwt_moles_biogas_produced(); //moles of biogas produced
      let moles_biogas_flared = moles_biogas*this.wwt_biog_fla/100; //moles of biogas flared
      let moles_ch4_flared    = moles_biogas_flared*this.wwt_ch4_biog/100; //moles of CH4 flared

      //combustion of 1 mol of CH4 produces 1 mol of CO2
      //CH4 + 2·O2 -> CO2 + 2·H2O
      //we do not account moles of CO2 already present into the biogas, because it is biogenic CO2
      let moles_co2_to_atmosphere = moles_ch4_flared; //moles of CO2
      let mass_co2_to_atmosphere = moles_co2_to_atmosphere*(44/1000); //kg of CO2

      let co2 = mass_co2_to_atmosphere; //kgCO2
      let n2o = 0;
      let ch4 = 0;
      let total = co2+ch4+n2o;
      return {total,co2,ch4,n2o};
    }

    //biogas valorized emissions
    wwt_KPI_GHG_biog_valorized(){
      let moles_biogas           = this.wwt_moles_biogas_produced(); //moles of biogas produced
      let moles_biogas_valorized = moles_biogas*this.wwt_biog_val/100; //moles of biogas valorized
      let moles_ch4_valorized    = moles_biogas_valorized*this.wwt_ch4_biog/100; //moles of CH4 valorized

      //combustion of 1 mol of CH4 produces 1 mol of CO2
      //CH4 + 2·O2 -> CO2 + 2·H2O
      //we do not account moles of CO2 already present into the biogas, because it is biogenic CO2
      let moles_co2_to_atmosphere = moles_ch4_valorized; //moles of CO2
      let mass_co2_to_atmosphere = moles_co2_to_atmosphere*(44/1000); //kg of CO2

      let co2 = mass_co2_to_atmosphere; //kgCO2
      let n2o = 0;
      let ch4 = 0;
      let total = co2+ch4+n2o;
      return {total,co2,ch4,n2o};
    }

    //biogas leaked emissions
    wwt_KPI_GHG_biog_leaked(){
      let moles_biogas = this.wwt_moles_biogas_produced(); //moles of biogas produced
      let moles_biogas_leaked = moles_biogas*this.wwt_biog_lkd/100; //moles of biogas leaked

      //we do not account moles of CO2 already present into the biogas, because it is biogenic CO2
      let moles_ch4_leaked = moles_biogas_leaked*this.wwt_ch4_biog/100; //moles of CH4 leaked
      let mass_ch4_to_atmosphere = moles_ch4_leaked*(16/1000); //kg of CH4 leaked

      let co2 = 0;
      let n2o = 0;
      let ch4 = mass_ch4_to_atmosphere*Cts.ct_ch4_eq.value; //kgCO2eq
      let total = co2+ch4+n2o;
      return {total,co2,ch4,n2o};
    }

    //ghg from sludge management
    wwt_KPI_GHG_slu(){
      let sources=[
        this.wwt_KPI_GHG_sludge_storage(),
        this.wwt_KPI_GHG_sludge_composting(),
        this.wwt_KPI_GHG_sludge_incineration(),
        this.wwt_KPI_GHG_sludge_land_application(),
        this.wwt_KPI_GHG_sludge_landfilling(),
        this.wwt_KPI_GHG_sludge_stockpilling(),
        this.wwt_KPI_GHG_sludge_transport(),
      ];

      //gases (numbers)
      let co2 = sources.map(s=>s.co2).sum();
      let ch4 = sources.map(s=>s.ch4).sum();
      let n2o = sources.map(s=>s.n2o).sum();

      //total
      let total = sources.map(s=>s.total).sum();
      return {total,co2,ch4,n2o};
    }

    //emissions from sludge storage
    wwt_KPI_GHG_sludge_storage(){
      let sludge_mass = this.wwt_mass_slu_sto; //kg of sludge
      let TVS         = this.wwt_slu_sto_TVS/100; //gTVS/gSludge
      let TVS_to_OC   = Cts.ct_VS_to_OC.value; //gOC/gTVS
      let OC_to_CH4   = Cts.ct_C_to_CH4_16_12.value; //gCH4/gOC
      let f_CH4       = this.wwt_slu_sto_f_CH4/100; //ratio for CH4 potential

      //max CH4 that could be released
      let ch4_potential = sludge_mass*TVS*TVS_to_OC*OC_to_CH4*f_CH4; //kgCH4 potential

      //emission factor
      let CH4_EF = this.wwt_slu_sto_EF/100; //gCH4 released / gCH4 potential

      //gases
      let co2   = 0;
      let n2o   = 0;
      let ch4   = ch4_potential*CH4_EF*Cts.ct_ch4_eq.value;
      let total = co2+ch4+n2o;
      return {total,co2,ch4,n2o};
    }

    //emissions from sludge composting
    wwt_KPI_GHG_sludge_composting(){
      let sludge_mass = this.wwt_mass_slu_comp; //kg of sludge
      let emissions_are_treated_or_piles_are_covered = this.wwt_slu_comp_emis_treated_or_piles_covered; //yes/no
      let solids_content_of_compost = this.wwt_slu_comp_solids_content; //%

      let TVS       = this.wwt_slu_comp_TVS/100; //gTVS/gSludge
      let N_cont    = this.wwt_slu_comp_N_cont/100; //gN/gSludge
      let TVS_to_OC = Cts.ct_VS_to_OC.value;  //gOC/gTVS
      let low_CN_EF = this.wwt_slu_comp_low_CN_EF; //0.015 kgN2O-N/kgN
      let up_EF     = this.wwt_slu_comp_uncovered_pile_EF; //0.025 kgCH4-C/kgC

      //gases
      let co2 = 0;
      let ch4 = (function(){
        if(emissions_are_treated_or_piles_are_covered){return 0}
        if(solids_content_of_compost>55){return 0}

        let OC_to_CH4 = Cts.ct_C_to_CH4_16_12.value; //1.33 gCH4/gOC
        return sludge_mass*TVS*TVS_to_OC*up_EF*OC_to_CH4*Cts.ct_ch4_eq.value;
      })();

      let n2o = (function(){
        let C_content = sludge_mass*TVS*TVS_to_OC; //kgC
        let N_content = sludge_mass*N_cont; //kgN
        let ratio_CN  = C_content/N_content||0;

        if(ratio_CN>30){return 0}
        if(solids_content_of_compost>55){return 0}

        return sludge_mass*N_cont*low_CN_EF*Cts.ct_N_to_N2O_44_28.value*Cts.ct_n2o_eq.value;
      })();

      let total = co2+ch4+n2o;
      return {total,co2,ch4,n2o};
    }

    //emissions from sludge incineration
    wwt_KPI_GHG_sludge_incineration(){
      let sludge_mass = this.wwt_mass_slu_inc;       //kg of sludge incinerated
      let Tf          = this.wwt_temp_inc;           //K
      let N_cont      = this.wwt_slu_inc_N_cont/100; //gN/gSludge
      let SNCR        = this.wwt_slu_inc_SNCR;       //yes/no

      //if Tf < 750ºC, use 750 ºC (1023 K)
      if(Tf < 1023){ Tf = 1023 }

      //gases
      let co2 = 0;
      let ch4 = (4.85e-5)*sludge_mass*Cts.ct_ch4_eq.value;
      let n2o = (function(){
        //n = % of total N that is emitted as N2O (suzuki et al 2003)
        let n = (161.3-0.14*Tf)/100; //gN2O/gN
        if(n<0) return 0;

        let emission = sludge_mass*N_cont*n*Cts.ct_n2o_eq.value; //kgCO2eq

        //increase N2O emissions by 20% if SNCR is used
        if(SNCR) emission *= 1.2;

        return emission;
      })();

      let total = co2+ch4+n2o;
      return {total,co2,ch4,n2o};
    }

    //emissions from sludge applied to land
    wwt_KPI_GHG_sludge_land_application(){
      let sludge_mass    = this.wwt_mass_slu_app; //kg sludge
      let solids_content = this.wwt_slu_la_solids_content; //%
      let TVS            = this.wwt_slu_la_TVS/100; //gTVS/gSludge
      let N_cont         = this.wwt_slu_la_N_cont/100; //gN/gSludge
      let TVS_to_OC      = Cts.ct_VS_to_OC.value; //gOC/gTVS
      let EF             = this.wwt_slu_la_EF; //gN2O-N/gN

      //gases
      let co2 = 0;
      let ch4 = 0;
      let n2o = (function(){
        //calculate ratio C:N
        let C_content = sludge_mass*TVS*TVS_to_OC; //kgC
        let N_content = sludge_mass*N_cont;        //kgN
        let ratio_CN  = C_content/N_content||0;

        if(ratio_CN>30){return 0}

        let emission = sludge_mass*N_cont*EF*Cts.ct_N_to_N2O_44_28.value*Cts.ct_n2o_eq.value;

        //if biosolids are >80%, N2O emissions are reduced by 50%
        if(solids_content>80) emission *= 0.5;

        return emission;
      })();

      let total = co2+ch4+n2o;
      return {total,co2,ch4,n2o};
    }

    //emissions from sludge used for landfilling
    wwt_KPI_GHG_sludge_landfilling(){
      let sludge_mass = this.wwt_mass_slu_land;

      let TVS         = this.wwt_slu_lf_TVS/100; //gTVS/gSludge
      let TVS_to_OC   = Cts.ct_VS_to_OC.value; //gOC/gTVS
      let uncertainty = this.wwt_slu_lf_uncertainty; //adimensional
      let OC_to_CH4   = Cts.ct_C_to_CH4_16_12.value; //kgCH4/kgC
      let CH4_in_gas  = this.wwt_slu_lf_CH4_in_gas/100; //%
      let DOCf        = this.wwt_slu_lf_DOCf/100; //%
      let decomp_3yr  = this.wwt_slu_lf_decomp_3yr/100; //%
      let MCF         = this.wwt_slu_lf_MCF; //methane correction for anaerobic managed landfills

      let N_cont      = this.wwt_slu_lf_N_cont/100; //gN/gSludge
      let low_CN_EF   = this.wwt_slu_lf_low_CN_EF; //0.015 kgN2O-N/kgN
      let N_to_N2O    = Cts.ct_N_to_N2O_44_28.value; //kgN2O/kgN2O-N

      let co2 = 0;
      let ch4 = (function(){
        return sludge_mass*TVS*TVS_to_OC*uncertainty*OC_to_CH4*CH4_in_gas*DOCf*decomp_3yr*MCF*Cts.ct_ch4_eq.value;
      })();
      let n2o = (function(){
        let C_cont = TVS*TVS_to_OC; //gOC/gSludge
        let ratio_CN = C_cont/N_cont||0; //gOC/gN
        if(ratio_CN>30){return 0;}
        return sludge_mass*N_cont*low_CN_EF*N_to_N2O*Cts.ct_n2o_eq.value;
      })();

      let total = co2+ch4+n2o;
      return {total,co2,ch4,n2o};
    }

    //emissions from sludge stockpiled
    wwt_KPI_GHG_sludge_stockpilling(){
      let sludge_mass = this.wwt_mass_slu_stock; //kg biosolids
      let sp_lifespan = this.wwt_slu_sp_lifespan; //years

      if(sp_lifespan<0) sp_lifespan=0;

      //integer part and decimal part (for example: 4.5 = 4 + 0.5)
      let sp_lifespan_int = Math.floor(sp_lifespan);
      let sp_lifespan_dec = sp_lifespan - sp_lifespan_int;

      //table of emission rates from Majumder et al., 2014 (table 3)
      let rates={
        //     <1 yo,  1-3 yo,   >3 yo, [kgCO2eq/kgSludge/year]
        ch4:[ 0.2e-3,  2.0e-3,  9.8e-3],
        n2o:[60.0e-3, 26.8e-3, 17.4e-3],
        co2:[30.1e-3, 30.5e-3, 10.1e-3],
      };

      //calculate emissions for 20 years then adapt to the real lifespan
      let emissions={ch4:[], n2o:[], co2:[]};

      //first year
      emissions.ch4[0] = sludge_mass*rates.ch4[0];
      emissions.n2o[0] = sludge_mass*rates.n2o[0];
      emissions.co2[0] = sludge_mass*rates.co2[0];

      //year 1 to 3
      for(let i=1;i<3;i++){
        emissions.ch4[i] = sludge_mass*rates.ch4[1];
        emissions.n2o[i] = sludge_mass*rates.n2o[1];
        emissions.co2[i] = sludge_mass*rates.co2[1];
      }

      //year 3 to 20
      for(let i=3;i<20;i++){
        emissions.ch4[i] = sludge_mass*rates.ch4[2];
        emissions.n2o[i] = sludge_mass*rates.n2o[2];
        emissions.co2[i] = sludge_mass*rates.co2[2];
      }

      //adapt emissions to real lifespan of stockpile (initialize to 0)
      let ch4 = 0;
      let n2o = 0;
      let co2 = 0;

      //sum whole years (integer part)
      for(let i=0; i < sp_lifespan_int; i++){
        ch4 += (emissions.ch4[i]||0);
        n2o += (emissions.n2o[i]||0);
        co2 += (emissions.co2[i]||0);
      }

      //sum decimal part
      ch4 += sp_lifespan_dec*(emissions.ch4[sp_lifespan_int]||0);
      n2o += sp_lifespan_dec*(emissions.n2o[sp_lifespan_int]||0);
      co2 += sp_lifespan_dec*(emissions.co2[sp_lifespan_int]||0);

      let total = co2+ch4+n2o;
      return {total,co2,ch4,n2o};
    }

    //emissions from sludge transport
    wwt_KPI_GHG_sludge_transport(){
      let vol   = this.wwt_vol_tslu;
      let fuel  = Tables.get_row('Fuel type',this.wwt_trck_typ);
      let co2   = vol*fuel.FD*fuel.NCV/1000*fuel.EFCO2;
      let n2o   = vol*fuel.FD*fuel.NCV/1000*fuel.EFN2O.vehicles*Cts.ct_n2o_eq.value;
      let ch4   = vol*fuel.FD*fuel.NCV/1000*fuel.EFCH4.vehicles*Cts.ct_ch4_eq.value;
      let total = co2+n2o+ch4;
      return {total,co2,n2o,ch4};
    }

    //emissions from water reuse transport
    wwt_KPI_GHG_reus_trck(){
      let vol   = this.wwt_reus_vol_trck;
      let fuel  = Tables.get_row('Fuel type',this.wwt_reus_trck_typ);
      let co2   = vol*fuel.FD*fuel.NCV/1000*fuel.EFCO2;
      let ch4   = vol*fuel.FD*fuel.NCV/1000*fuel.EFCH4.vehicles*Cts.ct_ch4_eq.value;
      let n2o   = vol*fuel.FD*fuel.NCV/1000*fuel.EFN2O.vehicles*Cts.ct_n2o_eq.value;
      let total = co2+n2o+ch4;
      return {total,co2,ch4,n2o};
    }

    //emissions from water discharged
    wwt_KPI_GHG_disc(){
      let co2   = 0;
      let ch4   = this.wwt_bod_effl*this.wwt_ch4_efac_dis*Cts.ct_ch4_eq.value;
      let n2o   = this.wwt_tn_effl *this.wwt_n2o_efac_dis*Cts.ct_N_to_N2O_44_28.value*Cts.ct_n2o_eq.value;
      let total = co2+ch4+n2o;
      return {total,co2,ch4,n2o};
    }

    /*
      Other methods
    */
    //convert volume of biogas produced to moles of biogas produced
    wwt_moles_biogas_produced(){
      //use PV=nRT formula
      //n = PV/RT
      //use normal conditions of pressure and temperature
      const P = 1.013e5; //Pa == N/m2 == J/m3
      let V = this.wwt_biog_pro; //m3
      const R = 8.31446261815324; //J/K·mol
      const T = 273.15; //K == 0ºC
      return P*V/(R*T); //"moles" of biogas produced
    }
};

/*
  TABLES
  Data structure for tabled values or dropdown menus
  used for two kinds of variables:
   1. nominal variables (strings) with magnitude==Option
   2. numeric variables inside Exceptions
*/
  let Tables={
    "Fuel type":[//            EF (kg/TJ)                                                                      FD (kg/L) NCV (TJ/Gg)
      {name:"Diesel",          EFCH4:{engines: 3,vehicles:3.9}, EFN2O:{engines:0.6,vehicles:3.9}, EFCO2:74100, FD:0.84,  NCV:43.0},
      {name:"Gasoline/Petrol", EFCH4:{engines: 3,vehicles:3.8}, EFN2O:{engines:0.6,vehicles:1.9}, EFCO2:69300, FD:0.74,  NCV:44.3},
      {name:"Natural Gas",     EFCH4:{engines:10,vehicles:92 }, EFN2O:{engines:0.1,vehicles:0.2}, EFCO2:56100, FD:0.75,  NCV:48.0},
    ],

    //ipcc 2019, table 6.3 (updated) EF (kgCH4/kgBOD)
    "type_of_water_body":[
      {name:"Water body undefined",                                                                   ch4_efac:0     },
      {name:"Discharge to aquatic environments (Tier 1)",                                             ch4_efac:0.068 },
      {name:"Discharge to aquatic environments other than reservoirs, lakes, and estuaries (Tier 2)", ch4_efac:0.021 },
      {name:"Discharge to reservoirs, lakes, and estuaries (Tier 2)",                                 ch4_efac:0.114 },
      {name:"Stagnant sewer or anaerobic water body",                                                 ch4_efac:0.3   },
      {name:"Flowing sewer (open or closed)",                                                         ch4_efac:0     },
      {name:"Soil infiltration",                                                                      ch4_efac:0     },
    ],

    "type_of_sewer":[
      {name:"Type of sewer undefined",                ch4_efac:0},
      {name:"Stagnant sewer or anaerobic water body", ch4_efac:0.3},
      {name:"Flowing sewer (open or closed)",         ch4_efac:0},
    ],

    //ipcc 2019, table 6.3 (updated) EF (kgCH4/kgBOD)
    "type_of_treatment":[
      {name:"Type of treatment undefined",                                  ch4_efac:0,     },
      {name:"Centralised, aerobic, treatment plant",                        ch4_efac:0.018, },
      {name:"Anaerobic Reactor - CH4 recovery not considered",              ch4_efac:0.48,  },
      {name:"Anaerobic Reactor - CH4 recovery considered",                  ch4_efac:0.14,  },
      {name:"Anaerobic shallow lagoon and facultative lagoons (<2m depth)", ch4_efac:0.12,  },
      {name:"Anaerobic deep lagoon (>2m depth)",                            ch4_efac:0.48,  },
      {name:"Anaerobic Lagoon covered",                                     ch4_efac:0,     },
      {name:"Wetlands - Surface flow",                                      ch4_efac:0.24,  },
      {name:"Wetlands - Horizontal subsurface flow",                        ch4_efac:0.06,  },
      {name:"Wetlands - Vertical subsurface flow",                          ch4_efac:0.006, },
      {name:"Activated Sludge - Well managed",                              ch4_efac:0,     },
      {name:"Activated Sludge - Minor poorly aerated zones",                ch4_efac:0.06,  },
      {name:"Activated Sludge - Some aerated zones",                        ch4_efac:0.12,  },
      {name:"Activated Sludge - Not well managed",                          ch4_efac:0.18,  },
      {name:"Aerated Lagoon",                                               ch4_efac:0.06,  },
      {name:"Trickling Filter",                                             ch4_efac:0.036, },
    ],

    "Type of onsite treatment":[
      {name:"Type of treatment undefined",                      ch4_efac:0.00,   bod_rmvd_as_sludge_estm:0.0,},
      {name:"Anaerobic Digester",                               ch4_efac:0.48,   bod_rmvd_as_sludge_estm:0.10,},
      {name:"Imhoff Tanks",                                     ch4_efac:0.48,   bod_rmvd_as_sludge_estm:0.10,},
      {name:"Anaerobic Reactors - CH4 recovery not considered", ch4_efac:0.48,   bod_rmvd_as_sludge_estm:0.10,},
      {name:"Anaerobic Reactors - CH4 recovery considered",     ch4_efac:0.14,   bod_rmvd_as_sludge_estm:0.10,},
      {name:"Stabilization Ponds (<2m depth)",                  ch4_efac:0.12,   bod_rmvd_as_sludge_estm:0.30,},
      {name:"Stabilization Ponds (>2m depth)",                  ch4_efac:0.48,   bod_rmvd_as_sludge_estm:0.10,},
      {name:"Sludge Drying Beds",                               ch4_efac:0.00,   bod_rmvd_as_sludge_estm:0.0,},
      {name:"Wetlands - surface flow",                          ch4_efac:0.24,   bod_rmvd_as_sludge_estm:0.30,},
      {name:"Wetlands - Horizontal subsurface flow",            ch4_efac:0.06,   bod_rmvd_as_sludge_estm:0.65,},
      {name:"Wetlands - Vertical subsurface flow",              ch4_efac:0.006,  bod_rmvd_as_sludge_estm:0.65,},
      {name:"Composting",                                       ch4_efac:0.0013, bod_rmvd_as_sludge_estm:0.0,},
      {name:"Activated Sludge (well managed)",                  ch4_efac:0.0000, bod_rmvd_as_sludge_estm:0.65,},
      {name:"Activated Sludge - minor poorly aerated zones",    ch4_efac:0.06,   bod_rmvd_as_sludge_estm:0.65,},
      {name:"Activated Sludge - Some aerated zones",            ch4_efac:0.12,   bod_rmvd_as_sludge_estm:0.65,},
      {name:"Activated Sludge - Not well managed",              ch4_efac:0.18,   bod_rmvd_as_sludge_estm:0.65,},
      {name:"Trickling Filter",                                 ch4_efac:0.036,  bod_rmvd_as_sludge_estm:0.65,},
    ],

    "N2O EF plants (Table 6.8A)":[
      {name:"Type of treatment undefined",           n2o_efac:0      },
      {name:"Centralised, aerobic, treatment plant", n2o_efac:0.016  },
      {name:"Anaerobic reactor",                     n2o_efac:0      },
      {name:"Anaerobic lagoons",                     n2o_efac:0      },
      {name:"Septic tank",                           n2o_efac:0      },
      {name:"Septic tank + land dispersal field",    n2o_efac:0.0045 },
      {name:"Latrine",                               n2o_efac:0      },
    ],

    "N2O EF effluent (Table 6.8A)":[
      {name:"Discharge undefined",                                                                              n2o_efac:0.000},
      {name:"Freshwater, estuarine, and marine discharge (Tier 1)",                                             n2o_efac:0.005},
      {name:"Nutrient-impacted and/or hypoxic freshwater, estuarine, and marine discharge (Tier 3, if needed)", n2o_efac:0.019},
    ],

    "type_of_treatment_KREM":[
      {name:"Mechanical treatment plants (primary sedimentation sludge)",                                                                                 K_rem:0.50},
      {name:"Aerobic treatment plants with primary treatment (mixed primary and secondary sludge, untreated or treated aerobically)",                     K_rem:0.80},
      {name:"Aerobic treatment plants with primary treatment and anaerobic sludge digestion (mixed primary and secondary sludge, treated anaerobically)", K_rem:1.00},
      {name:"Aerobic wastewater treatment plants without separate primary treatment",                                                                     K_rem:1.16},
    ],

    "WW treatment organics removal fractions (centralised) (Table 6.6B and 6.10C)":[
      {name:"Untreated systems",                                                     bod_effl:1,    N_effl:1.00},
      {name:"Primary (mechanical treatment plants)",                                 bod_effl:0.60, N_effl:0.90},
      {name:"Primary + Secondary (biological treatment plants)",                     bod_effl:0.15, N_effl:0.60},
      {name:"Primary + Secondary + Tertiary (advanced biological treatment plants)", bod_effl:0.10, N_effl:0.20},
    ],

    "WW treatment organics removal fractions (onsite) (Table 6.6B and 6.10C)":[
      {name:"Untreated systems",                                                                        bod_rmvd:0,     N_effl:1.00 },
      {name:"Septic tank/septic system",                                                                bod_rmvd:0.625, N_effl:0.85 },
      {name:"Septic tank/septic system + land dispersal field",                                         bod_rmvd:0.625, N_effl:0.32 },
      {name:"Latrines – Dry climate, groundwater table lower than latrine, small family (3–5 persons)", bod_rmvd:0.1,   N_effl:0.88 },
      {name:"Latrines – Dry climate, groundwater table lower than latrine, communal (many users)",      bod_rmvd:0.5,   N_effl:0.88 },
      {name:"Latrines – Wet climate/flush water use, groundwater table higher than latrine",            bod_rmvd:0.7,   N_effl:0.88 },
    ],

    //Andreoli et al table 2.2
    "Sludge characteristics in each stage of the treatment process":[
      {name:"Type of treatment undefined",                                     gSS_inh_day:0},
      {name:"Primary treatment (conventional)",                                gSS_inh_day:(35+45)/2},
      {name:"Primary treatment (septic tanks)",                                gSS_inh_day:(20+30)/2},
      {name:"Facultative pond",                                                gSS_inh_day:(20+25)/2},
      {name:"Anaerobic pond – facultative pond (anaerobic pond)",              gSS_inh_day:(20+45)/2},
      {name:"Anaerobic pond – facultative pond (facultative pond)",            gSS_inh_day:( 6+10)/2},
      {name:"Anaerobic pond – facultative pond (total)",                       gSS_inh_day:(26+55)/2},
      {name:"Facultative aerated lagoon",                                      gSS_inh_day:( 8+13)/2},
      {name:"Complete-mix aerat.lagoon – sedim. pond",                         gSS_inh_day:(11+13)/2},
      {name:"Septic tank+anaerobic filter (septic tank)",                      gSS_inh_day:(20+30)/2},
      {name:"Septic tank+anaerobic filter (anaerobic filter)",                 gSS_inh_day:( 7+9 )/2},
      {name:"Septic tank+anaerobic filter (total)",                            gSS_inh_day:(27+39)/2},
      {name:"Conventional activated sludge (primary sludge)",                  gSS_inh_day:(35+45)/2},
      {name:"Conventional activated sludge (secondary sludge)",                gSS_inh_day:(25+35)/2},
      {name:"Conventional activated sludge (mixed sludge)",                    gSS_inh_day:(60+80)/2},
      {name:"Activated sludge extended aeration",                              gSS_inh_day:(40+45)/2},
      {name:"High rate trickling filter (primary sludge)",                     gSS_inh_day:(35+45)/2},
      {name:"High rate trickling filter (secondary sludge)",                   gSS_inh_day:(20+30)/2},
      {name:"High rate trickling filter (mixed sludge)",                       gSS_inh_day:(55+75)/2},
      {name:"Submerged aerated biofilter (primary sludge)",                    gSS_inh_day:(35+45)/2},
      {name:"Submerged aerated biofilter (secondary sludge)",                  gSS_inh_day:(25+35)/2},
      {name:"Submerged aerated biofilter (mixed sludge)",                      gSS_inh_day:(60+80)/2},
      {name:"UASB Reactor",                                                    gSS_inh_day:(12+18)/2},
      {name:"UASB+activated sludge (anaerobic sludge (UASB))",                 gSS_inh_day:(12+18)/2},
      {name:"UASB+activated sludge (aerobic sludge (activated sludge))",       gSS_inh_day:( 8+14)/2},
      {name:"UASB+activated sludge (mixed sludge)",                            gSS_inh_day:(20+32)/2},
      {name:"UASB+aerobic biofilm reactor (anaerobic sludge (UASB))",          gSS_inh_day:(12+18)/2},
      {name:"UASB+aerobic biofilm reactor (aerobic sludge (aerobic reactor))", gSS_inh_day:( 6+12)/2},
      {name:"UASB+aerobic biofilm reactor (mixed sludge)",                     gSS_inh_day:(18+30)/2},
    ],

    "Type of sludge disposed":[
      {name:"Type of sludge disposed", f_ch4:0,  N_cont:0, TVS:0 },
      {name:"Non-digested",            f_ch4:53, N_cont:3, TVS:70},
      {name:"Digested",                f_ch4:6,  N_cont:4, TVS:51},
    ],

    //for land application and landfilling
    "Type of faecal sludge":[
      {name:"Type of faecal sludge undefined", N_content:0.00, TVS:0,  total_solids:0.00},
      {name:"Untreated faecal sludge",         N_content:0.24, TVS:70, total_solids:0.04},
      {name:"Treated faecal sludge",           N_content:3.00, TVS:40, total_solids:0.22},
      {name:"Pit humus",                       N_content:4.00, TVS:65, total_solids:0.07},
      {name:"Dehydrated faeces",               N_content:3.00, TVS:70, total_solids:0.27},
      {name:"Compost",                         N_content:3.00, TVS:80, total_solids:0.08},
      {name:"Septic tank sludge",              N_content:0.03, TVS:60, total_solids:0.02},
    ],

    "Type of landfill":[
      {name:"Landfill",                     MCF:1},
      {name:"Landfill (with gas recovery)", MCF:0.02},
      {name:"Landfill (flaring)",           MCF:0},
    ],

    //f_la: gN transformed to gN2O
    "Soil type":[
      {name:"Soil type undefined",         f_la:0.000},
      {name:"Fine-Textured (>30% clay)",   f_la:0.023},
      {name:"Coarse-Textured (<30% clay)", f_la:0.005},
    ],

    "Type of containment":[
      {name:"Containment undefined",                                          ch4_efac:0,      ch4_efac_flooding:0,     BOD_conc_FS:0,    fs_density:0   },
      {name:"No containment (open defecation)",                               ch4_efac:0.027,  ch4_efac_flooding:0.027, BOD_conc_FS:67.8, fs_density:1400},
      {name:"Pit latrine without flush water (lined or unlined) – household", ch4_efac:0.06,   ch4_efac_flooding:0.42,  BOD_conc_FS:67.8, fs_density:1400},
      {name:"Pit latrine without flush water (lined or unlined) – communal",  ch4_efac:0.3,    ch4_efac_flooding:0.42,  BOD_conc_FS:67.8, fs_density:1400},
      {name:"Pit latrine with flush water use (lined or unlined)",            ch4_efac:0.42,   ch4_efac_flooding:0.42,  BOD_conc_FS:67.8, fs_density:1400},
      {name:"Septic tank (with or without dispersal field)",                  ch4_efac:0.3,    ch4_efac_flooding:0.42,  BOD_conc_FS:1.35, fs_density:1100},
      {name:"Fully lined tank without flush water use – not water tight",     ch4_efac:0.3,    ch4_efac_flooding:0.42,  BOD_conc_FS:67.8, fs_density:1400},
      {name:"Fully lined tank without flush water use – water tight",         ch4_efac:0.42,   ch4_efac_flooding:0.42,  BOD_conc_FS:67.8, fs_density:1400},
      {name:"Fully lined tank with flush water use - water tight or untight", ch4_efac:0.42,   ch4_efac_flooding:0.42,  BOD_conc_FS:67.8, fs_density:1400},
      {name:"Urine Diverting Dry Toilet (UDDT)",                              ch4_efac:0.0,    ch4_efac_flooding:0.42,  BOD_conc_FS:67.8, fs_density:1400},
      {name:"Composting Toilet",                                              ch4_efac:0.0013, ch4_efac_flooding:0.42,  BOD_conc_FS:67.8, fs_density:1400},
      {name:"Imhoff tank",                                                    ch4_efac:0.48,   ch4_efac_flooding:0.42,  BOD_conc_FS:67.8, fs_density:1400},
    ],

    "Yes/No":[
      {name:"no"},
      {name:"yes"},
    ],
  };

  //get row object by "table" (string) and "index" (integer)
  Tables.get_row=function(table, index){
    let t=Tables[table]; //array

    //check if table exists
    if(!t       ) throw `Table "${table}" does not exist`;
    if(!t[index]) throw `Table.${table}[${index}] does not exist`;

    //checks passed: return row
    return t[index];
  };

  //table titles and/or descriptions
  Tables.get_table_description=function(table_name){
    return {
      "type_of_water_body":    "EFCH4 for Type of Water Body (Table 6.3)",
      "type_of_sewer":         "EFCH4 for Type of Sewer (Table 6.3)",
      "type_of_treatment":     "EFCH4 for Type of Treatment (Table 6.3)",
      "type_of_treatment_KREM":"Removal of organic component from wastewater as sludge (KREM) according to treatment type (Table 6.6A)",
    }[table_name]||table_name;
  };

/*
  CONSTANTS
*/
let Cts={
  ct_ch4_eq:{
    value:34,
    descr:"Conversion of CH4 emissions to CO2 equivalent emissions",
    unit:"kgCO2eq/kgCH4",
  },
  ct_n2o_eq:{
    value:298,
    descr:"Conversion of N2O emissions to CO2 equivalent emissions",
    unit:"kgCO2eq/kgN2O",
  },
  ct_C_to_CH4_16_12:{
    value:16/12,
    descr:"Organic C to CH4 conversion factor (=16/12)",
    unit:"gCH4/gOC"
  },
  ct_N_to_N2O_44_28:{
    value:44/28,
    descr:"N2O-N to N2O conversion factor (=44/28)",
    unit:"gN2O/gN2O-N",
  },
  ct_C_to_CO2_44_12:{
    value:44/12,
    descr:"C to CO2 conversion (=44/12)",
    unit:"gCO2/gC",
  },
  ct_VS_to_OC:{
    value:0.56,
    descr:"Organic Carbon content in Volatile Solids",
    unit:"gOC/gVS",
  },
};

/*
  Assessment example:
    - create an assessment
    - add 1 industry
    - modify some inputs
    - calculate total emissions
*/
let a = new Assessment();
a.add_industry();
a.industries[0].wwt_nrg_cons = 500;
a.industries[0].wwt_conv_kwh = 0.5;
console.log(a.TotalGHG());
