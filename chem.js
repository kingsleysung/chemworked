// ChemSteps shared chemistry logic: atomic masses, formula parsing, significant figures.
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

  // Escapes text before it goes into innerHTML (error messages can echo what the user typed).
  function esc(s) {
    return String(s).replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
  }

  const Chem = {
    esc, ELEMENTS, parseFormula, molarMass,
    parseNumberString, countSigFigs, lastDecimalPlace, formatSigFigs, roundToPlace, sigFigArithmetic,
    toSuperscript, plain,
  };

  if (typeof module !== "undefined" && module.exports) module.exports = Chem;
  else root.Chem = Chem;
})(typeof window !== "undefined" ? window : globalThis);
