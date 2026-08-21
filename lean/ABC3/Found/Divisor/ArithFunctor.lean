/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.ArithOrd
import Mathlib.RingTheory.Ideal.GoingUp
import Mathlib.NumberTheory.RamificationInertia.Valuation
import Mathlib.NumberTheory.RamificationInertia.Inertia

/-!
# `Φ` と `B` の数体についての関手性(鎖 `arith` の `arith-functor`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.113。

原文 (FrdI p.113):
> v = ord(Fv), all but a finite number of which are zero]. Note, moreover, that

## ★★何を作るか

原文は「`Φ`, `B`, および準同型 `B → Φ^gp` は数体 `F` について関手的である」とだけ書く。
その中身は、有限拡大 `L ⊆ M` に対する

  `Φ(L)^gp → Φ(M)^gp`

である。★これは `ord(F_v) → ord(F'_w)`(`F_v ⊆ F'_w` が誘導する写像)であり、
**対数の模型では「局所次数倍」というただ 1 つの規則**になる:

| 素点 | 係数 | 拡大での倍率 |
|---|---|---|
| 非アルキメデス | `ord_v(x)·log(N v)` | `e·f` |
| アルキメデス | `-mult(v)·log(v x)` | `mult(W)/mult(v)` |

★どちらも**局所次数 `[M_W : L_v]`** である。★★したがって

  `arithExtend d := (V ↦ localDeg(V) · d(res V))`

の 1 本で済む。

## ★在庫(2026-08-20 に測った)

* `IsDedekindDomain.HeightOneSpectrum.valuation_liesOver` —— `v.valuation x ^ e = w.valuation (algebraMap x)`
* `Ideal.absNorm_eq_pow_inertiaDeg_of_liesOver` —— `N W = (N v)^f`
* `NumberField.InfinitePlace.comap_apply` / `mult_comap_le`
* `IsDedekindDomain.primesOver_finite` —— ファイバーの有限性
* `NumberField.FinitePlace.maximalIdeal_injective`
-/

namespace ABC3.Found.Divisor

open NumberField IsDedekindDomain Finsupp

universe u

variable {L M : Type u} [Field L] [Field M] [NumberField L] [NumberField M] [Algebra L M]

/-! ## ★1. 素点の制限 -/

/-- ★高さ 1 の素点の制限(`𝓞 M` の素イデアルを `𝓞 L` へ落とす)。 -/
noncomputable def resHOS (V : HeightOneSpectrum (𝓞 M)) : HeightOneSpectrum (𝓞 L) where
  asIdeal := V.asIdeal.under (𝓞 L)
  isPrime := Ideal.IsPrime.under _ _
  ne_bot := Ideal.IsIntegral.comap_ne_bot (𝓞 L) V.ne_bot

instance liesOver_resHOS (V : HeightOneSpectrum (𝓞 M)) :
    V.asIdeal.LiesOver (resHOS (L := L) V).asIdeal := ⟨rfl⟩

/-- ★有限素点の制限。 -/
noncomputable def resFin (W : FinitePlace M) : FinitePlace L :=
  FinitePlace.mk (resHOS (L := L) W.maximalIdeal)

/-- ★無限素点の制限。 -/
noncomputable def resInf (W : InfinitePlace M) : InfinitePlace L :=
  W.comap (algebraMap L M)

/-- ★★**素点の制限** `V(M) → V(L)`。 -/
noncomputable def resPlace : ArithPlace M → ArithPlace L :=
  Sum.elim (fun W => Sum.inl (resFin (L := L) W)) (fun W => Sum.inr (resInf (L := L) W))

@[simp] theorem resPlace_inl (W : FinitePlace M) :
    resPlace (L := L) (Sum.inl W) = Sum.inl (resFin (L := L) W) := rfl

@[simp] theorem resPlace_inr (W : InfinitePlace M) :
    resPlace (L := L) (Sum.inr W) = Sum.inr (resInf (L := L) W) := rfl

theorem maximalIdeal_resFin (W : FinitePlace M) :
    (resFin (L := L) W).maximalIdeal = resHOS (L := L) W.maximalIdeal :=
  FinitePlace.maximalIdeal_mk _

theorem logAbsNorm_resFin (W : FinitePlace M) :
    logAbsNorm (resFin (L := L) W)
      = Real.log ((Ideal.absNorm (resHOS (L := L) W.maximalIdeal).asIdeal : ℕ) : ℝ) := by
  rw [logAbsNorm, maximalIdeal_resFin]

/-! ## ★2. 局所次数 -/

/-- ★分岐指数 `e(W|v)`。 -/
noncomputable def ramIdx (W : FinitePlace M) : ℕ :=
  (resHOS (L := L) W.maximalIdeal).asIdeal.ramificationIdx W.maximalIdeal.asIdeal

/-- ★剰余次数 `f(W|v)`。 -/
noncomputable def inertDeg (W : FinitePlace M) : ℕ :=
  (resHOS (L := L) W.maximalIdeal).asIdeal.inertiaDeg W.maximalIdeal.asIdeal

theorem ramIdx_pos (W : FinitePlace M) : 0 < ramIdx (L := L) W :=
  Nat.pos_of_ne_zero (Ideal.IsDedekindDomain.ramificationIdx_ne_zero_of_liesOver
    W.maximalIdeal.asIdeal (resHOS (L := L) W.maximalIdeal).ne_bot)

theorem inertDeg_pos (W : FinitePlace M) : 0 < inertDeg (L := L) W :=
  Ideal.inertiaDeg_pos _ _

theorem mult_eq_one_or_two {K : Type*} [Field K] [NumberField K] (w : InfinitePlace K) :
    w.mult = 1 ∨ w.mult = 2 := by
  unfold InfinitePlace.mult
  split <;> simp

theorem mult_resInf_dvd (W : InfinitePlace M) :
    (resInf (L := L) W).mult ∣ W.mult := by
  have hle : (resInf (L := L) W).mult ≤ W.mult :=
    InfinitePlace.mult_comap_le (algebraMap L M) W
  rcases mult_eq_one_or_two (resInf (L := L) W) with h1 | h1
  · rw [h1]; exact one_dvd _
  · rcases mult_eq_one_or_two W with h2 | h2
    · rw [h1, h2] at hle; omega
    · rw [h1, h2]

/-- ★★**局所次数 `[M_W : L_v]`** —— 非アルキメデスで `e·f`、アルキメデスで `mult(W)/mult(v)`。 -/
noncomputable def localDeg : ArithPlace M → ℕ :=
  Sum.elim (fun W => ramIdx (L := L) W * inertDeg (L := L) W)
    (fun W => W.mult / (resInf (L := L) W).mult)

@[simp] theorem localDeg_inl (W : FinitePlace M) :
    localDeg (L := L) (Sum.inl W) = ramIdx (L := L) W * inertDeg (L := L) W := rfl

@[simp] theorem localDeg_inr (W : InfinitePlace M) :
    localDeg (L := L) (Sum.inr W) = W.mult / (resInf (L := L) W).mult := rfl

theorem localDeg_pos (V : ArithPlace M) : 0 < localDeg (L := L) V := by
  cases V with
  | inl W => exact Nat.mul_pos (ramIdx_pos W) (inertDeg_pos W)
  | inr W =>
      refine Nat.div_pos ?_ ?_
      · exact InfinitePlace.mult_comap_le (algebraMap L M) W
      · rcases mult_eq_one_or_two (resInf (L := L) W) with h | h <;> omega

theorem localDeg_mul_mult (W : InfinitePlace M) :
    localDeg (L := L) (Sum.inr W) * (resInf (L := L) W).mult = W.mult := by
  rw [localDeg_inr]
  exact Nat.div_mul_cancel (mult_resInf_dvd W)

/-! ## ★3. 局所での恒等式 -/

/-- ★★**`ord_W(x) = e(W|v)·ord_v(x)`**(`x ∈ L`)。

★`HeightOneSpectrum.valuation_liesOver` から。値群 `Multiplicative ℤ` の中で
`a^e = b` を `toAdd` で読むだけである。 -/
theorem ordFin_extend (W : FinitePlace M) {x : L} (hx : x ≠ 0) :
    ordFin W.maximalIdeal (algebraMap L M x)
      = (ramIdx (L := L) W : ℤ) * ordFin (resHOS (L := L) W.maximalIdeal) x := by
  set V := W.maximalIdeal with hV
  set w := resHOS (L := L) V with hw
  set e := ramIdx (L := L) W with he
  have hval : (w.valuation L x) ^ e = V.valuation M (algebraMap L M x) :=
    HeightOneSpectrum.valuation_liesOver M w V x
  have hwx : w.valuation L x ≠ 0 := (Valuation.ne_zero_iff _).mpr hx
  have hVx : V.valuation M (algebraMap L M x) ≠ 0 := by
    rw [← hval]
    exact pow_ne_zero _ hwx
  have hcoe : ((WithZero.unzero hwx) ^ e : Multiplicative ℤ) = WithZero.unzero hVx := by
    refine WithZero.coe_inj.mp ?_
    rw [WithZero.coe_unzero hVx, ← hval]
    push_cast
    rw [WithZero.coe_unzero hwx]
  rw [ordFin, ordFin, dif_neg hVx, dif_neg hwx, ← hcoe]
  show -((WithZero.unzero hwx) ^ e).toAdd = (e : ℤ) * -(WithZero.unzero hwx).toAdd
  have hpow : ((WithZero.unzero hwx) ^ e).toAdd = (e : ℤ) * (WithZero.unzero hwx).toAdd := by
    simp
  rw [hpow]
  ring

/-- ★★**`log(N W) = f(W|v)·log(N v)`**。 -/
theorem logAbsNorm_extend (W : FinitePlace M) :
    logAbsNorm W = (inertDeg (L := L) W : ℝ) * logAbsNorm (resFin (L := L) W) := by
  have h : Ideal.absNorm W.maximalIdeal.asIdeal
      = Ideal.absNorm (resHOS (L := L) W.maximalIdeal).asIdeal ^ (inertDeg (L := L) W) :=
    Ideal.absNorm_eq_pow_inertiaDeg_of_liesOver _ _
      (resHOS (L := L) W.maximalIdeal).isPrime (resHOS (L := L) W.maximalIdeal).ne_bot
  rw [logAbsNorm, logAbsNorm_resFin, h]
  push_cast
  exact Real.log_pow _ _

/-- ★★★★**局所の恒等式** —— `arithPlaceLog` は拡大で**局所次数倍**になる。 -/
theorem arithPlaceLog_extend {x : L} (hx : x ≠ 0) (V : ArithPlace M) :
    arithPlaceLog (algebraMap L M x) V
      = (localDeg (L := L) V : ℝ) * arithPlaceLog x (resPlace (L := L) V) := by
  have hax : algebraMap L M x ≠ 0 := (map_ne_zero_iff _ (algebraMap L M).injective).mpr hx
  cases V with
  | inl W =>
      rw [resPlace_inl, localDeg_inl]
      have h1 : arithPlaceLog (algebraMap L M x) (Sum.inl W)
          = (ordFin W.maximalIdeal (algebraMap L M x) : ℝ) * logAbsNorm W := by
        have h := arithPlaceLog_finite W.maximalIdeal (algebraMap L M x) hax
        rw [FinitePlace.mk_maximalIdeal] at h
        exact h
      have h2 : arithPlaceLog x (Sum.inl (resFin (L := L) W))
          = (ordFin (resHOS (L := L) W.maximalIdeal) x : ℝ)
            * logAbsNorm (resFin (L := L) W) := by
        rw [logAbsNorm_resFin]
        exact arithPlaceLog_finite (resHOS (L := L) W.maximalIdeal) x hx
      rw [h1, h2, ordFin_extend W hx, logAbsNorm_extend (L := L) W]
      push_cast
      ring
  | inr W =>
      rw [resPlace_inr]
      have hW : W (algebraMap L M x) = (resInf (L := L) W) x :=
        (InfinitePlace.comap_apply W (algebraMap L M) x).symm
      show -(W.mult : ℝ) * Real.log (W (algebraMap L M x))
        = (localDeg (L := L) (Sum.inr W) : ℝ)
          * (-((resInf (L := L) W).mult : ℝ) * Real.log ((resInf (L := L) W) x))
      rw [hW]
      have hm : ((localDeg (L := L) (Sum.inr W) : ℕ) : ℝ) * ((resInf (L := L) W).mult : ℝ)
          = (W.mult : ℝ) := by
        exact_mod_cast congrArg (fun n : ℕ => (n : ℝ)) (localDeg_mul_mult (L := L) W)
      rw [← hm]
      ring

/-! ## ★4. ファイバーの有限性 -/

theorem finite_fiber_fin (v : FinitePlace L) :
    {W : FinitePlace M | resFin (L := L) W = v}.Finite := by
  classical
  haveI : v.maximalIdeal.asIdeal.IsMaximal :=
    Ideal.IsPrime.isMaximal v.maximalIdeal.isPrime v.maximalIdeal.ne_bot
  have hinj : Function.Injective (fun W : FinitePlace M => W.maximalIdeal.asIdeal) :=
    HeightOneSpectrum.asIdeal_injective.comp FinitePlace.maximalIdeal_injective
  refine Set.Finite.of_finite_image ?_ (Function.Injective.injOn hinj)
  refine Set.Finite.subset
    (IsDedekindDomain.primesOver_finite v.maximalIdeal.asIdeal (𝓞 M)) ?_
  rintro _ ⟨W, hW, rfl⟩
  refine ⟨W.maximalIdeal.isPrime, ⟨?_⟩⟩
  have h1 : resHOS (L := L) W.maximalIdeal = v.maximalIdeal := by
    have h := congrArg FinitePlace.maximalIdeal hW
    rwa [maximalIdeal_resFin] at h
  exact (congrArg HeightOneSpectrum.asIdeal h1).symm

theorem finite_fiber (v : ArithPlace L) :
    {V : ArithPlace M | resPlace (L := L) V = v}.Finite := by
  classical
  cases v with
  | inl w =>
      have hsub : {V : ArithPlace M | resPlace (L := L) V = Sum.inl w}
          ⊆ Sum.inl '' {W : FinitePlace M | resFin (L := L) W = w} := by
        rintro V hV
        cases V with
        | inl W =>
            refine ⟨W, ?_, rfl⟩
            simpa using hV
        | inr W => simp at hV
      exact Set.Finite.subset ((finite_fiber_fin w).image Sum.inl) hsub
  | inr w =>
      refine Set.Finite.subset (Set.finite_range (Sum.inr : InfinitePlace M → ArithPlace M)) ?_
      rintro V hV
      cases V with
      | inl W => simp at hV
      | inr W => exact ⟨W, rfl⟩

theorem finite_preimage (T : Finset (ArithPlace L)) :
    (resPlace (L := L) (M := M) ⁻¹' (T : Set (ArithPlace L))).Finite := by
  classical
  have h : resPlace (L := L) (M := M) ⁻¹' (T : Set (ArithPlace L))
      = ⋃ v ∈ (T : Set (ArithPlace L)), {V : ArithPlace M | resPlace (L := L) V = v} := by
    ext V
    simp [Set.mem_preimage]
  rw [h]
  exact T.finite_toSet.biUnion (fun v _ => finite_fiber v)

/-! ## ★5. 拡大写像 -/

theorem arithExtend_support_aux (d : ArithPlace L →₀ ℝ) :
    ∀ V : ArithPlace M, (localDeg (L := L) V : ℝ) * d (resPlace (L := L) V) ≠ 0 →
      V ∈ (finite_preimage (L := L) (M := M) d.support).toFinset := by
  intro V hV
  refine (Set.Finite.mem_toFinset _).mpr (Set.mem_preimage.mpr ?_)
  refine Finset.mem_coe.mpr (Finsupp.mem_support_iff.mpr ?_)
  intro h0
  rw [h0, mul_zero] at hV
  exact hV rfl

/-- ★★★★**`Φ(L)^gp → Φ(M)^gp`** —— 局所次数倍。 -/
noncomputable def arithExtend :
    (ArithPlace L →₀ ℝ) →+ (ArithPlace M →₀ ℝ) where
  toFun d := Finsupp.onFinset (finite_preimage (L := L) (M := M) d.support).toFinset
    (fun V => (localDeg (L := L) V : ℝ) * d (resPlace (L := L) V))
    (arithExtend_support_aux d)
  map_zero' := by
    ext V
    show (localDeg (L := L) V : ℝ) * (0 : ArithPlace L →₀ ℝ) (resPlace (L := L) V) = 0
    simp
  map_add' d e := by
    ext V
    show (localDeg (L := L) V : ℝ) * (d + e) (resPlace (L := L) V)
      = (localDeg (L := L) V : ℝ) * d (resPlace (L := L) V)
        + (localDeg (L := L) V : ℝ) * e (resPlace (L := L) V)
    rw [Finsupp.add_apply]
    ring

@[simp] theorem arithExtend_apply (d : ArithPlace L →₀ ℝ) (V : ArithPlace M) :
    arithExtend (L := L) d V = (localDeg (L := L) V : ℝ) * d (resPlace (L := L) V) := rfl

/-- ★★**`arithExtend` は `Φ^gp` を `Φ^gp` へ移す**。

★非アルキメデス成分は `(e·f·n)·log(N v) = (e·n)·log(N W)` になる。 -/
theorem arithExtend_mem_arithDivGroup {d : ArithPlace L →₀ ℝ} (hd : d ∈ arithDivGroup L) :
    arithExtend (L := L) d ∈ arithDivGroup M := by
  intro W
  obtain ⟨n, hn⟩ := hd (resFin (L := L) W)
  refine ⟨(ramIdx (L := L) W : ℤ) * n, ?_⟩
  show (localDeg (L := L) (Sum.inl W) : ℝ) * d (resPlace (L := L) (Sum.inl W))
      = ((((ramIdx (L := L) W : ℤ) * n : ℤ)) : ℝ) * logAbsNorm W
  rw [resPlace_inl]
  rw [show d (Sum.inl (resFin (L := L) W))
      = (n : ℝ) * logAbsNorm (resFin (L := L) W) from hn]
  rw [localDeg_inl, logAbsNorm_extend (L := L) W]
  push_cast
  ring

/-- ★★**`arithExtend` は有効因子を有効因子へ移す**。 -/
theorem arithExtend_mem_arithEff {d : ArithPlace L →₀ ℝ} (hd : d ∈ arithEff L) :
    arithExtend (L := L) d ∈ arithEff M := by
  refine ⟨arithExtend_mem_arithDivGroup hd.1, fun V => ?_⟩
  show (0 : ℝ) ≤ (localDeg (L := L) V : ℝ) * d (resPlace (L := L) V)
  exact mul_nonneg (by positivity) (hd.2 _)

/-- ★★★★★**`B → Φ^gp` の関手性** ——
`div^M(x) = arithExtend (div^L(x))`(`x ∈ L^×`)。

原文 (FrdI p.113):
> F. Thus, if F is a [not necessarily finite] Galois extension of F, G -/
theorem arithExtend_arithDiv (x : Lˣ) :
    arithDiv (L := M) (Units.map (algebraMap L M).toMonoidHom x)
      = arithExtend (L := L) (arithDiv x) := by
  ext V
  have hx : (x : L) ≠ 0 := x.ne_zero
  show arithPlaceLog ((algebraMap L M) (x : L)) V
    = (localDeg (L := L) V : ℝ) * arithPlaceLog (x : L) (resPlace (L := L) V)
  exact arithPlaceLog_extend hx V

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Example 6.3` の「`Φ`, `B`, `B → Φ^gp` は数体について関手的」。 -/
def arithExtend.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — Φ・B・B → Φ^gp が数体 F について関手的",
    sectionId := "frdi-example-6-3" }

/-! ## ★★同型に沿った局所次数は 1 —— non-dilating の入力(2026-08-21)

`Theorem 6.4, (i)` の non-dilating は
`isNonDilating_of_primary_sharp`(`MonoidTransport.lean`)の入力 `hfix`
「自己同型は素点を素点へ**係数 1 で**移す」に落ちる。
★`arithExtend d V = localDeg V · d (resPlace V)` なので、要るのは
**同型に沿った `localDeg V = 1`** ただ 1 つである。

★★**finrank 経由は塞がっている** —— `letI := σ.toAlgebra` を置いても
`Module.finrank L L` は標準の module 構造に解決されるので、
`sum_ramification_inertia` に載せる筋は書けない。
★そこで **`e` と `f` を直接**計算する:

| 量 | 根拠 |
|---|---|
| `e = 1` | `map σ (comap σ P) = P`(全射)＋ Dedekind 環で `P² < P` |
| `f = 1` | 剰余体の間の射が全射(商どうしなので型が違い、instance の衝突が無い) |
| 無限素点 | `InfinitePlace.isReal_comap_iff` で `mult` が保たれる |
-/

attribute [local instance] Ideal.Quotient.field

omit [NumberField L] in
/-- ★★**全射なら分岐指数は 1**。

★`P` は Dedekind 環の 0 でない素イデアルなので `P² < P`、
したがって `P ≤ P²` は起こらない。 -/
theorem ramIdx_eq_one_of_surjective
    (hs : Function.Surjective (algebraMap (𝓞 L) (𝓞 M))) (W : FinitePlace M) :
    ramIdx (L := L) W = 1 := by
  haveI : (W.maximalIdeal.asIdeal).IsPrime := W.maximalIdeal.isPrime
  have hne : W.maximalIdeal.asIdeal ≠ ⊥ := W.maximalIdeal.ne_bot
  have hmap : Ideal.map (algebraMap (𝓞 L) (𝓞 M))
      (resHOS (L := L) W.maximalIdeal).asIdeal = W.maximalIdeal.asIdeal :=
    Ideal.map_comap_of_surjective _ hs _
  refine Ideal.ramificationIdx_spec ?_ ?_
  · rw [hmap, pow_one]
  · rw [hmap]
    intro h
    have hlt := Ideal.pow_succ_lt_pow (P := W.maximalIdeal.asIdeal) hne 1
    rw [pow_one] at hlt
    exact absurd (le_antisymm hlt.le h) (ne_of_lt hlt)

/-- ★★**全射なら剰余次数は 1**。

★剰余体の間の射も全射で、**商どうしなので型が違い instance の衝突が無い**。 -/
theorem inertDeg_eq_one_of_surjective
    (hs : Function.Surjective (algebraMap (𝓞 L) (𝓞 M))) (W : FinitePlace M) :
    inertDeg (L := L) W = 1 := by
  haveI hpm : (resHOS (L := L) W.maximalIdeal).asIdeal.IsMaximal :=
    Ideal.IsPrime.isMaximal (resHOS (L := L) W.maximalIdeal).isPrime
      (resHOS (L := L) W.maximalIdeal).ne_bot
  rw [inertDeg, Ideal.inertiaDeg_algebraMap]
  refine (finrank_eq_one_iff_of_nonzero' (1 : _) one_ne_zero).mpr ?_
  intro w
  obtain ⟨y, rfl⟩ := Ideal.Quotient.mk_surjective w
  obtain ⟨x, rfl⟩ := hs y
  refine ⟨Ideal.Quotient.mk _ x, ?_⟩
  rw [Algebra.smul_def, mul_one]
  rfl

omit [NumberField L] [NumberField M] in
/-- ★★**同型に沿えば無限素点の `mult` は保たれる**(実素点は実素点へ)。 -/
theorem mult_resInf_of_bijective (hb : Function.Bijective (algebraMap L M))
    (W : InfinitePlace M) : (resInf (L := L) W).mult = W.mult := by
  classical
  have hiff : (resInf (L := L) W).IsReal ↔ W.IsReal :=
    InfinitePlace.isReal_comap_iff (f := RingEquiv.ofBijective (algebraMap L M) hb) (w := W)
  unfold InfinitePlace.mult
  by_cases hr : W.IsReal
  · simp [hiff.mpr hr, hr]
  · have h2 : ¬ (resInf (L := L) W).IsReal := fun h => hr (hiff.mp h)
    simp [h2, hr]

/-- ★★★★**同型に沿った局所次数は 1**。

★これが `Theorem 6.4, (i)` の non-dilating の入力 `hfix` の中身である。 -/
theorem localDeg_eq_one_of_bijective (hb : Function.Bijective (algebraMap L M))
    (hs : Function.Surjective (algebraMap (𝓞 L) (𝓞 M))) (V : ArithPlace M) :
    localDeg (L := L) V = 1 := by
  cases V with
  | inl W => rw [localDeg_inl, ramIdx_eq_one_of_surjective hs W,
      inertDeg_eq_one_of_surjective hs W, one_mul]
  | inr W =>
      rw [localDeg_inr, mult_resInf_of_bijective hb W,
        Nat.div_self InfinitePlace.mult_pos]

/-- ★★★locator —— 同型に沿った局所次数が 1 であること
(`Theorem 6.4, (i)` の non-dilating の入力)。 -/
def localDeg_eq_one_of_bijective.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — 自己同型に沿った局所次数は 1(non-dilating の入力)",
    sectionId := "frdi-example-6-3" }

end ABC3.Found.Divisor
