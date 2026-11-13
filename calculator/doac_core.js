/* DOAC dose by Cockcroft-Gault vs eGFR — core logic (shared by web calculator and tests).
 * Independent JavaScript re-implementation of egfr.py + doac_rules.py.
 * Creatinine in µmol/L; 1 mg/dL = 88.4 µmol/L. */
(function (root) {
  "use strict";
  var U = 88.4;

  function cg(age, wt, scr, female) {            // mL/min, actual body weight
    var v = (140 - age) * wt / (72 * (scr / U));
    return female ? v * 0.85 : v;
  }
  function ckdepi2009(age, scr, female) {
    var s = scr / U, k = female ? 0.7 : 0.9, a = female ? -0.329 : -0.411, r = s / k;
    var v = 141 * Math.pow(Math.min(r, 1), a) * Math.pow(Math.max(r, 1), -1.209) * Math.pow(0.993, age);
    return female ? v * 1.018 : v;
  }
  function ckdepi2021(age, scr, female) {
    var s = scr / U, k = female ? 0.7 : 0.9, a = female ? -0.241 : -0.302, r = s / k;
    var v = 142 * Math.pow(Math.min(r, 1), a) * Math.pow(Math.max(r, 1), -1.200) * Math.pow(0.9938, age);
    return female ? v * 1.012 : v;
  }
  function ekfcQ(age, female) {                   // Pottel 2021 Fig. 1: formula 2–25 y, fixed >25 y
    if (age <= 25) {
      var l = female
        ? 3.080 + 0.177 * age - 0.223 * Math.log(age) - 0.00596 * age * age + 0.0000686 * age * age * age
        : 3.200 + 0.259 * age - 0.543 * Math.log(age) - 0.00763 * age * age + 0.0000790 * age * age * age;
      return Math.exp(l);
    }
    return female ? 62 : 80;
  }
  function ekfc(age, scr, female) {
    var r = scr / ekfcQ(age, female);
    var v = 107.3 * Math.pow(r, r < 1 ? -0.322 : -1.132);
    return age > 40 ? v * Math.pow(0.990, age - 40) : v;
  }
  function bsaDuBois(wt, ht) { return 0.007184 * Math.pow(wt, 0.425) * Math.pow(ht, 0.725); }
  function bsaMosteller(wt, ht) { return Math.sqrt(wt * ht / 3600); }
  function deindex(v, bsa) { return v * bsa / 1.73; }

  /* Dose levels: 0 not recommended / CI, 1 reduced, 1.5 dabigatran individual choice,
     2 standard, 3 (US edoxaban) not used because CrCl > 95. */
  var EU = {
    rivaroxaban: function (rf) { return rf < 15 ? 0 : rf < 50 ? 1 : 2; },
    edoxaban: function (rf, p) {
      var l = rf < 15 ? 0 : rf <= 50 ? 1 : 2;
      return (l === 2 && (p.weight <= 60 || p.pgp)) ? 1 : l;
    },
    dabigatran: function (rf, p) {
      if (rf < 30) return 0;
      if (p.age >= 80 || p.verapamil) return 1;
      if (rf <= 50 || (p.age >= 75 && p.age < 80)) return 1.5;
      return 2;
    },
    apixaban: function (rf, p) {
      if (rf < 15) return 0;
      if (rf < 30) return 1;
      var n = (p.age >= 80 ? 1 : 0) + (p.weight <= 60 ? 1 : 0) + (p.scr >= 133 ? 1 : 0);
      return n >= 2 ? 1 : 2;
    }
  };
  var US = {
    rivaroxaban: function (rf) { return rf <= 50 ? 1 : 2; },
    edoxaban: function (rf) { return rf < 15 ? 0 : rf <= 50 ? 1 : rf <= 95 ? 2 : 3; },
    dabigatran: function (rf) { return rf < 15 ? 0 : rf <= 30 ? 1 : 2; },
    apixaban: function (rf, p) {
      var n = (p.age >= 80 ? 1 : 0) + (p.weight <= 60 ? 1 : 0) + (p.scr >= 133 ? 1 : 0);
      return n >= 2 ? 1 : 2;
    }
  };
  var LABELS = {
    EU: {
      rivaroxaban: { 0: "Not recommended", 1: "15 mg once daily", 2: "20 mg once daily" },
      edoxaban: { 0: "Not recommended", 1: "30 mg once daily", 2: "60 mg once daily" },
      dabigatran: { 0: "Contraindicated", 1: "110 mg twice daily", 1.5: "150 or 110 mg twice daily (individual)", 2: "150 mg twice daily" },
      apixaban: { 0: "Not recommended", 1: "2.5 mg twice daily", 2: "5 mg twice daily" }
    },
    US: {
      rivaroxaban: { 1: "15 mg once daily", 2: "20 mg once daily" },
      edoxaban: { 0: "Not recommended", 1: "30 mg once daily", 2: "60 mg once daily", 3: "Do not use (CrCl > 95)" },
      dabigatran: { 0: "No recommendation", 1: "75 mg twice daily", 2: "150 mg twice daily" },
      apixaban: { 1: "2.5 mg twice daily", 2: "5 mg twice daily" }
    }
  };

  function estimates(p) {
    var f = !!p.female, bD = bsaDuBois(p.weight, p.height), bM = bsaMosteller(p.weight, p.height);
    var e = { CG: cg(p.age, p.weight, p.scr, f), BSA_DuBois: bD, BSA_Mosteller: bM };
    [["CKDEPI09", ckdepi2009], ["CKDEPI21", ckdepi2021], ["EKFC", ekfc]].forEach(function (x) {
      var v = x[1](p.age, p.scr, f);
      e[x[0]] = v; e[x[0] + "_abs"] = deindex(v, bD); e[x[0] + "_absM"] = deindex(v, bM);
    });
    return e;
  }
  function doses(p, region) {
    var rules = region === "US" ? US : EU, e = estimates(p), out = {};
    Object.keys(rules).forEach(function (drug) {
      out[drug] = {};
      Object.keys(e).forEach(function (k) {
        if (k.indexOf("BSA") === 0) return;
        out[drug][k] = rules[drug](e[k], p);
      });
    });
    return { estimates: e, doses: out };
  }

  var api = { cg: cg, ckdepi2009: ckdepi2009, ckdepi2021: ckdepi2021, ekfc: ekfc, ekfcQ: ekfcQ,
              bsaDuBois: bsaDuBois, bsaMosteller: bsaMosteller, deindex: deindex,
              EU: EU, US: US, LABELS: LABELS, estimates: estimates, doses: doses, VERSION: "1.0.0" };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.DOAC = api;
})(this);
