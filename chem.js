// ChemWorked shared chemistry logic: atomic masses, formula parsing, significant figures.
// Works in the browser (window.Chem) and in Node (module.exports) so it can be tested.

(function (root) {
  // Standard atomic weights (IUPAC abridged values; mass number of the most stable isotope for radioactive elements).
  const ELEMENT_DATA =
    "H,Hydrogen,1.008|He,Helium,4.0026|Li,Lithium,6.94|Be,Beryllium,9.0122|B,Boron,10.81|C,Carbon,12.011|" +
    "N,Nitrogen,14.007|O,Oxygen,15.999|F,Fluorine,18.998|Ne,Neon,20.180|Na,Sodium,22.990|Mg,Magnesium,24.305|" +
    "Al,Aluminium,26.982|Si,Silicon,28.085|P,Phosphorus,30.974|S,Sulfur,32.06|Cl,Chlorine,35.45|Ar,Argon,39.948|" +
    "K,Potassium,39.098|Ca,Calcium,40.078|Sc,Scandium,44.956|Ti,Titanium,47.867|V,Vanadium,50.942|Cr,Chromium,51.996|" +
    "Mn,Manganese,54.938|Fe,Iron,55.845|Co,Cobalt,58.933|Ni,Nickel,58.693|Cu,Copper,63.546|Zn,Zinc,65.38|" +
    "Ga,Gallium,69.723|Ge,Germanium,72.630|As,Arsenic,74.922|Se,Selenium,78.971|Br,Bromine,79.904|Kr,Krypton,83.798|" +
    "Rb,Rubidium,85.468|Sr,Strontium,87.62|Y,Yttrium,88.906|Zr,Zirconium,91.224|Nb,Niobium,92.906|Mo,Molybdenum,95.95|" +
    "Tc,Technetium,98|Ru,Ruthenium,101.07|Rh,Rhodium,102.91|Pd,Palladium,106.42|Ag,Silver,107.87|Cd,Cadmium,112.41|" +
    "In,Indium,114.82|Sn,Tin,118.71|Sb,Antimony,121.76|Te,Tellurium,127.60|I,Iodine,126.90|Xe,Xenon,131.29|" +
    "Cs,Caesium,132.91|Ba,Barium,137.33|La,Lanthanum,138.91|Ce,Cerium,140.12|Pr,Praseodymium,140.91|Nd,Neodymium,144.24|" +
    "Pm,Promethium,145|Sm,Samarium,150.36|Eu,Europium,151.96|Gd,Gadolinium,157.25|Tb,Terbium,158.93|Dy,Dysprosium,162.50|" +
    "Ho,Holmium,164.93|Er,Erbium,167.26|Tm,Thulium,168.93|Yb,Ytterbium,173.05|Lu,Lutetium,174.97|Hf,Hafnium,178.49|" +
    "Ta,Tantalum,180.95|W,Tungsten,183.84|Re,Rhenium,186.21|Os,Osmium,190.23|Ir,Iridium,192.22|Pt,Platinum,195.08|" +
    "Au,Gold,196.97|Hg,Mercury,200.59|Tl,Thallium,204.38|Pb,Lead,207.2|Bi,Bismuth,208.98|Po,Polonium,209|" +
    "At,Astatine,210|Rn,Radon,222|Fr,Francium,223|Ra,Radium,226|Ac,Actinium,227|Th,Thorium,232.04|" +
    "Pa,Protactinium,231.04|U,Uranium,238.03|Np,Neptunium,237|Pu,Plutonium,244";

  const ELEMENTS = {};
  ELEMENT_DATA.split("|").forEach((row) => {
    const [symbol, name, mass] = row.split(",");
    ELEMENTS[symbol] = { name, mass: parseFloat(mass) };
  });

  // ---------- Formula parsing ----------

  const OPEN = { "(": ")", "[": "]", "{": "}" };
  const CLOSE = new Set([")", "]", "}"]);

  // Parses formulas like H2O, Ca(OH)2, K4[Fe(CN)6], CuSO4·5H2O (also CuSO4*5H2O or CuSO4.5H2O).
  // Returns { counts: {symbol: n}, order: [symbols in order of first appearance] }.
  function parseFormula(input) {
    const formula = input.replace(/\s+/g, "");
    if (!formula) throw new Error("Type a chemical formula, like H2O or Ca(OH)2.");

    const counts = {};
    const order = [];
    const add = (sym, n) => {
      if (!(sym in counts)) { counts[sym] = 0; order.push(sym); }
      counts[sym] += n;
    };

    // Hydrates and adducts: split on a dot, then each part may start with a coefficient (5H2O).
    const parts = formula.split(/[·•*.]/);
    parts.forEach((part) => {
      if (!part) throw new Error("There is an empty part around the dot (·). Example of a hydrate: CuSO4·5H2O");
      const m = part.match(/^(\d+)(.*)$/);
      const coef = m ? parseInt(m[1], 10) : 1;
      const body = m ? m[2] : part;
      if (!body) throw new Error(`"${part}" has a number but no formula after it.`);
      const { group, index } = parseGroup(body, 0, null);
      if (index !== body.length) throw new Error(`Unexpected "${body[index]}" — check your brackets.`);
      group.order.forEach((sym) => add(sym, group.counts[sym] * coef));
    });

    return { counts, order };
  }

  function readNumber(s, i) {
    let j = i;
    while (j < s.length && /\d/.test(s[j])) j++;
    return { n: j > i ? parseInt(s.slice(i, j), 10) : 1, index: j };
  }

  function parseGroup(s, i, closer) {
    const counts = {};
    const order = [];
    const add = (sym, n) => {
      if (!(sym in counts)) { counts[sym] = 0; order.push(sym); }
      counts[sym] += n;
    };

    while (i < s.length) {
      const c = s[i];
      if (c in OPEN) {
        const inner = parseGroup(s, i + 1, OPEN[c]);
        const num = readNumber(s, inner.index);
        inner.group.order.forEach((sym) => add(sym, inner.group.counts[sym] * num.n));
        i = num.index;
      } else if (CLOSE.has(c)) {
        if (c !== closer) throw new Error(`The bracket "${c}" doesn't match an opening bracket.`);
        if (!order.length) throw new Error("There are empty brackets in the formula.");
        return { group: { counts, order }, index: i + 1 };
      } else if (/[A-Z]/.test(c)) {
        let sym = c;
        if (i + 1 < s.length && /[a-z]/.test(s[i + 1])) sym += s[i + 1];
        if (!ELEMENTS[sym]) {
          if (sym.length === 2 && ELEMENTS[c]) {
            throw new Error(`"${sym}" isn't an element. Element symbols are a capital letter plus at most one lowercase letter.`);
          }
          throw new Error(`"${sym}" isn't an element symbol.`);
        }
        const num = readNumber(s, i + sym.length);
        if (num.n === 0) throw new Error(`"${sym}0" — a subscript can't be zero.`);
        add(sym, num.n);
        i = num.index;
      } else if (/[a-z]/.test(c)) {
        throw new Error(`Element symbols start with a capital letter. Did you mean "${c.toUpperCase()}${s.slice(i + 1, i + 2)}"? (For example CO is carbon monoxide, Co is cobalt.)`);
      } else if (/\d/.test(c)) {
        throw new Error(`A number appears before any element ("${c}"). Subscripts go after the element, like H2O.`);
      } else {
        throw new Error(`"${c}" isn't allowed in a formula.`);
      }
    }

    if (closer) throw new Error(`A bracket was opened but never closed — add "${closer}".`);
    return { group: { counts, order }, index: i };
  }

  // Returns rows for a step-by-step table plus the total molar mass.
  function molarMass(input) {
    const { counts, order } = parseFormula(input);
    let total = 0;
    const rows = order.map((sym) => {
      const el = ELEMENTS[sym];
      const subtotal = el.mass * counts[sym];
      total += subtotal;
      return { symbol: sym, name: el.name, count: counts[sym], mass: el.mass, subtotal };
    });
    rows.forEach((r) => { r.percent = (r.subtotal / total) * 100; });
    return { rows, total };
  }

  // ---------- Significant figures ----------

  // Accepts "0.00450", "1200", "1.20e3", "1.20E-3", "1.20x10^3", "1.20×10^3", "-3.0".
  function parseNumberString(raw) {
    const s = String(raw).trim().replace(/\s+/g, "").replace(/,/g, "");
    let m = s.match(/^([+-]?)(\d*\.?\d*)(?:[eE]([+-]?\d+)|[x×*]10\^?\(?([+\-−]?\d+)\)?)?$/);
    if (!m || !/\d/.test(m[2])) return null;
    const exp = m[3] !== undefined ? parseInt(m[3], 10)
      : m[4] !== undefined ? parseInt(m[4].replace("−", "-"), 10) : 0;
    const mantissa = m[2];
    return {
      sign: m[1] === "-" ? -1 : 1,
      mantissa,
      exponent: exp,
      scientific: m[3] !== undefined || m[4] !== undefined,
      value: (m[1] === "-" ? -1 : 1) * parseFloat(mantissa) * Math.pow(10, exp),
    };
  }

  // Counts significant figures and marks which digits are significant.
  // Returns { count, ambiguous, digits: [{char, significant}] } or null if invalid.
  function countSigFigs(raw) {
    const p = parseNumberString(raw);
    if (!p) return null;
    const chars = p.mantissa.split("");
    const hasDecimal = p.mantissa.includes(".");
    const digitIdx = chars.map((c, i) => (/\d/.test(c) ? i : -1)).filter((i) => i >= 0);
    const firstNonZero = digitIdx.find((i) => chars[i] !== "0");

    if (firstNonZero === undefined) {
      // The number is zero: only zeros after the decimal point are measured.
      const after = hasDecimal ? p.mantissa.split(".")[1].length : 0;
      const count = Math.max(after, 1);
      return {
        count,
        ambiguous: false,
        isZero: true,
        digits: chars.map((c) => ({ char: c, significant: false })),
        parsed: p,
      };
    }

    let last;
    let ambiguous = false;
    if (hasDecimal) {
      last = digitIdx[digitIdx.length - 1];
    } else {
      last = [...digitIdx].reverse().find((i) => chars[i] !== "0");
      ambiguous = last !== digitIdx[digitIdx.length - 1];
    }

    let count = 0;
    const digits = chars.map((c, i) => {
      const significant = /\d/.test(c) && i >= firstNonZero && i <= last;
      if (significant) count++;
      return { char: c, significant };
    });
    return { count, ambiguous, isZero: false, digits, parsed: p };
  }

  // Decimal place of the last significant digit: 2 means hundredths, -2 means hundreds.
  function lastDecimalPlace(raw) {
    const info = countSigFigs(raw);
    if (!info) return null;
    const p = info.parsed;
    const [intPart, fracPart = ""] = p.mantissa.split(".");
    let place;
    if (p.mantissa.includes(".")) {
      place = fracPart.length;
    } else {
      const trailing = intPart.match(/0*$/)[0].length;
      place = info.isZero ? 0 : -trailing;
    }
    return place - p.exponent;
  }

  // Formats a number with exactly n significant figures, using ×10ⁿ when trailing zeros would be ambiguous.
  function formatSigFigs(x, n) {
    if (!isFinite(x)) return "—";
    if (n < 1) n = 1;
    if (x === 0) return n > 1 ? "0." + "0".repeat(n - 1) : "0";
    const nudged = x + Math.sign(x) * Math.abs(x) * 1e-12; // round halves away from zero (school convention)
    let s = nudged.toPrecision(n);
    if (/e/.test(s)) {
      const [mant, e] = s.split("e");
      const exp = parseInt(e, 10);
      // Small exponents read better as plain decimals (0.00450 instead of 4.50 × 10⁻³) unless very small.
      if (exp < 0 && exp >= -4) return Number(s).toFixed(n - 1 - exp);
      return `${mant} × 10${toSuperscript(exp)}`;
    }
    if (!s.includes(".")) {
      const digits = s.replace("-", "");
      if (/0$/.test(digits) && digits.length === n) s += "."; // "100." shows the zeros are significant
    }
    return s;
  }

  function toSuperscript(n) {
    const map = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };
    return String(n).split("").map((c) => map[c] || c).join("");
  }

  // Rounds x so its last digit sits at the given decimal place, returning a formatted string.
  function roundToPlace(x, place) {
    if (x === 0) return place > 0 ? (0).toFixed(place) : "0";
    const nudged = x + Math.sign(x) * Math.abs(x) * 1e-12;
    const scale = Math.pow(10, place);
    const rounded = (Math.sign(nudged) * Math.round(Math.abs(nudged) * scale)) / scale;
    if (rounded === 0) return place > 0 ? (0).toFixed(place) : "0";
    // Rounding can carry into a new digit (9.96 → 10.0), so count digits after rounding.
    const exponent = Math.floor(Math.log10(Math.abs(rounded)) + 1e-12);
    return formatSigFigs(rounded, exponent + place + 1);
  }

  // Sig-fig arithmetic for two measured values.
  function sigFigArithmetic(aRaw, op, bRaw) {
    const a = countSigFigs(aRaw);
    const b = countSigFigs(bRaw);
    if (!a || !b) return null;
    const av = a.parsed.value;
    const bv = b.parsed.value;
    let exact;
    switch (op) {
      case "+": exact = av + bv; break;
      case "-": exact = av - bv; break;
      case "×": exact = av * bv; break;
      case "÷":
        if (bv === 0) throw new Error("You can't divide by zero.");
        exact = av / bv; break;
      default: return null;
    }
    if (op === "+" || op === "-") {
      const pa = lastDecimalPlace(aRaw);
      const pb = lastDecimalPlace(bRaw);
      const place = Math.min(pa, pb);
      return { exact, answer: roundToPlace(exact, place), rule: "decimal", placeA: pa, placeB: pb, place, a, b };
    }
    const sf = Math.min(a.count, b.count);
    return { exact, answer: formatSigFigs(exact, sf), rule: "sigfigs", sf, a, b };
  }

  // Plain readable number for intermediate steps (not rounded for sig figs).
  function plain(x, maxDecimals = 6) {
    if (!isFinite(x)) return "—";
    if (x !== 0 && (Math.abs(x) < 1e-4 || Math.abs(x) >= 1e9)) {
      const [m, e] = x.toExponential(5).split("e");
      return `${parseFloat(m)} × 10${toSuperscript(parseInt(e, 10))}`;
    }
    return String(parseFloat(x.toFixed(maxDecimals)));
  }

  // ---------- Equation balancing ----------

  // Exact fractions with BigInt so balancing never suffers from rounding.
  const babs = (a) => (a < 0n ? -a : a);
  const bgcd = (a, b) => { a = babs(a); b = babs(b); while (b) [a, b] = [b, a % b]; return a; };
  function frac(n, d = 1n) {
    if (d < 0n) { n = -n; d = -d; }
    const g = bgcd(n, d) || 1n;
    return { n: n / g, d: d / g };
  }
  const fsub = (a, b) => frac(a.n * b.d - b.n * a.d, a.d * b.d);
  const fmul = (a, b) => frac(a.n * b.n, a.d * b.d);
  const fdiv = (a, b) => frac(a.n * b.d, a.d * b.n);

  // Balances equations like "Fe + O2 -> Fe2O3". Returns { reactants, products, elements, matrix }.
  function balanceEquation(input) {
    const text = input.replace(/\s+/g, " ").trim();
    const sides = text.split(/\s*(?:<?-+>|→|⟶|⇌|<=>|=)\s*/);
    if (sides.length !== 2 || !sides[0] || !sides[1]) {
      throw new Error('Write the equation with an arrow between reactants and products, like "H2 + O2 -> H2O".');
    }
    const parseSide = (side) => side.split(/\s*\+\s*/).map((raw) => {
      const formula = raw.trim().replace(/^\d+\s*/, ""); // drop any coefficient the user typed
      if (!formula) throw new Error("There's an empty spot next to a + sign.");
      return { formula, ...parseFormula(formula) };
    });
    const reactants = parseSide(sides[0]);
    const products = parseSide(sides[1]);
    const species = [...reactants, ...products];

    const elements = [];
    species.forEach((s) => s.order.forEach((el) => { if (!elements.includes(el)) elements.push(el); }));
    elements.forEach((el) => {
      const left = reactants.some((s) => s.counts[el]);
      const right = products.some((s) => s.counts[el]);
      if (!left || !right) throw new Error(`${el} appears on only one side, so this equation can't be balanced. Check your formulas.`);
    });

    // One row per element: reactant counts positive, product counts negative. Solve matrix · x = 0.
    const rows = elements.map((el) => species.map((s, j) =>
      frac(BigInt((s.counts[el] || 0) * (j < reactants.length ? 1 : -1)))));
    const n = species.length;
    const pivots = [];
    let r = 0;
    for (let c = 0; c < n && r < rows.length; c++) {
      let p = rows.findIndex((row, i) => i >= r && row[c].n !== 0n);
      if (p === -1) continue;
      [rows[r], rows[p]] = [rows[p], rows[r]];
      const pv = rows[r][c];
      rows[r] = rows[r].map((v) => fdiv(v, pv));
      for (let i = 0; i < rows.length; i++) {
        if (i !== r && rows[i][c].n !== 0n) {
          const f = rows[i][c];
          rows[i] = rows[i].map((v, k) => fsub(v, fmul(f, rows[r][k])));
        }
      }
      pivots.push(c);
      r++;
    }
    const free = [...Array(n).keys()].filter((c) => !pivots.includes(c));
    if (free.length === 0) throw new Error("This equation can't be balanced. Check that the formulas are right.");
    if (free.length > 1) throw new Error("This equation can be balanced in more than one way (it's really two reactions combined). Try splitting it up.");

    const f = free[0];
    const x = Array(n).fill(null);
    x[f] = frac(1n);
    pivots.forEach((c, i) => { x[c] = frac(-rows[i][f].n, rows[i][f].d); });
    const lcm = x.reduce((acc, v) => (acc * v.d) / bgcd(acc, v.d), 1n);
    let ints = x.map((v) => (v.n * lcm) / v.d);
    const g = ints.reduce((acc, v) => bgcd(acc, v), 0n);
    ints = ints.map((v) => v / g);
    if (ints.every((v) => v <= 0n)) ints = ints.map((v) => -v);
    if (ints.some((v) => v <= 0n)) throw new Error("This equation can't be balanced with these formulas. Check for a typo or a missing product.");

    const coefs = ints.map(Number);
    reactants.forEach((s, i) => { s.coef = coefs[i]; });
    products.forEach((s, i) => { s.coef = coefs[reactants.length + i]; });
    const tally = elements.map((el) => ({
      element: el,
      left: reactants.reduce((t, s) => t + s.coef * (s.counts[el] || 0), 0),
      right: products.reduce((t, s) => t + s.coef * (s.counts[el] || 0), 0),
    }));
    return { reactants, products, elements, tally };
  }

  // ---------- Empirical formula ----------

  // items: [{ symbol, mass }] where mass is grams or percent (percent = grams in a 100 g sample).
  function empiricalFormula(items) {
    if (items.length < 2) throw new Error("Enter at least two elements.");
    const rows = items.map(({ symbol, mass }) => {
      const el = ELEMENTS[symbol];
      if (!el) throw new Error(`"${symbol}" isn't an element symbol. Use a capital first letter, like Cl.`);
      if (!(mass > 0)) throw new Error(`Enter a positive amount for ${symbol}.`);
      return { symbol, mass, atomicMass: el.mass, moles: mass / el.mass };
    });
    const smallest = Math.min(...rows.map((r) => r.moles));
    rows.forEach((r) => { r.ratio = r.moles / smallest; });
    let multiplier = 1;
    for (let m = 1; m <= 8; m++) {
      if (rows.every((r) => Math.abs(r.ratio * m - Math.round(r.ratio * m)) < 0.1)) { multiplier = m; break; }
      if (m === 8) multiplier = null;
    }
    if (multiplier === null) throw new Error("The mole ratios don't come out close to whole numbers. Check your data.");
    rows.forEach((r) => { r.subscript = Math.round(r.ratio * multiplier); });
    const formula = rows.map((r) => r.symbol + (r.subscript > 1 ? r.subscript : "")).join("");
    const empiricalMass = rows.reduce((t, r) => t + r.subscript * r.atomicMass, 0);
    return { rows, smallest, multiplier, formula, empiricalMass };
  }

  // Escapes text before it goes into innerHTML (error messages can echo what the user typed).
  // Two-tailed Student's t at 95% confidence, indexed by degrees of freedom.
  const T95 = [NaN, 12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228, 2.201, 2.179, 2.160, 2.145,
    2.131, 2.120, 2.110, 2.101, 2.093, 2.086, 2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042];
  function t95(df) {
    if (df < 1) return NaN;
    if (df <= 30) return T95[df];
    if (df <= 40) return 2.021;
    if (df <= 60) return 2.000;
    if (df <= 120) return 1.980;
    return 1.960;
  }

  // Least-squares calibration line y = m·x + b, with the standard uncertainty of an unknown read off the line
  // (the usual analytical chemistry formula, e.g. Harris "Quantitative Chemical Analysis", ch. 4).
  // points: [{x, y}], unknownYs: replicate signals of the unknown (may be empty).
  function calibration(points, unknownYs = []) {
    const n = points.length;
    if (n < 3) return { error: "Enter at least 3 standards to fit a line and estimate its uncertainty." };
    const xbar = points.reduce((s, p) => s + p.x, 0) / n;
    const ybar = points.reduce((s, p) => s + p.y, 0) / n;
    let Sxx = 0, Sxy = 0, Syy = 0, sumX2 = 0;
    for (const p of points) {
      Sxx += (p.x - xbar) ** 2;
      Sxy += (p.x - xbar) * (p.y - ybar);
      Syy += (p.y - ybar) ** 2;
      sumX2 += p.x ** 2;
    }
    if (Sxx === 0) return { error: "All the standards have the same concentration, so a line can't be fit." };
    const m = Sxy / Sxx;
    const b = ybar - m * xbar;
    const rows = points.map((p) => ({ ...p, dx: p.x - xbar, dy: p.y - ybar, resid: p.y - (m * p.x + b) }));
    const SSres = rows.reduce((s, r) => s + r.resid ** 2, 0);
    const sy = Math.sqrt(SSres / (n - 2));
    const sm = sy / Math.sqrt(Sxx);
    const sb = sy * Math.sqrt(sumX2 / (n * Sxx));
    const r2 = Syy === 0 ? NaN : 1 - SSres / Syy;
    const fit = { n, xbar, ybar, Sxx, Sxy, Syy, sumX2, m, b, rows, SSres, sy, sm, sb, r2, df: n - 2, t: t95(n - 2) };

    const k = unknownYs.length;
    if (k > 0 && m !== 0) {
      const yu = unknownYs.reduce((s, y) => s + y, 0) / k;
      const x = (yu - b) / m;
      const sx = (sy / Math.abs(m)) * Math.sqrt(1 / k + 1 / n + (yu - ybar) ** 2 / (m ** 2 * Sxx));
      const xs = points.map((p) => p.x);
      fit.unknown = { k, yu, x, sx, ci: fit.t * sx, outside: x < Math.min(...xs) || x > Math.max(...xs) };
    }
    return fit;
  }

  // Round a value and its uncertainty together: uncertainty to 2 sig figs, value to the same decimal place.
  function withUncertainty(x, u) {
    if (!isFinite(x) || !isFinite(u) || u <= 0) return { value: plain(x), unc: plain(u) };
    const place = Math.floor(Math.log10(u)) - 1; // decimal place of the 2nd sig fig of u
    const fix = (v) => (place < 0 ? (Math.round(v / 10 ** place) * 10 ** place).toFixed(-place) : String(Math.round(v / 10 ** place) * 10 ** place));
    return { value: fix(x), unc: fix(u) };
  }

  function esc(s) {
    return String(s).replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
  }

  const Chem = {
    esc, ELEMENTS, parseFormula, molarMass, balanceEquation, empiricalFormula,
    parseNumberString, countSigFigs, lastDecimalPlace, formatSigFigs, roundToPlace, sigFigArithmetic,
    toSuperscript, plain, t95, calibration, withUncertainty,
  };

  if (typeof module !== "undefined" && module.exports) module.exports = Chem;
  else root.Chem = Chem;
})(typeof window !== "undefined" ? window : globalThis);
