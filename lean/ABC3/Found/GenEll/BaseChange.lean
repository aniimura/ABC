import ABC3.Found.GenEll.ArithDiv
import ABC3.Found.GenEll.FinitePlaceRel
import ABC3.Found.GenEll.InfinitePlaceRel
import Mathlib.NumberTheory.RamificationInertia.Inertia

/-!
# [GenEll] §1 —— 正規化次数の**底変換不変性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

## ★★これが何を塞ぐか —— `check.mjs` 冒頭 B5 の穴

`Skeleton/GenEll/Heights.lean` の `degNormalized_base_change` は
**`Interface` の `PulledBackClassData.base_change_invariant` を取り出しているだけ**で、
中身を証明していない(`.needs` に自分でそう書いてある)。

★本ファイルは**その中身を `ADiv` の層で実際に証明する**。

> **`deg_K(a_K) = [K:F]·deg_F(a)`**、したがって **`deg_K(a_K) = deg_F(a)`**(正規化次数)

★★これは `Interface` の仮説を消すものではない(あちらは `X` 上の直線束の話で、
こちらは `Spec(𝓞_F)` 上の算術因子の話である)。
**しかし `degNormalized` が「正規化する理由」として掲げていた不変性は、
ここで初めて仮説でなく定理になった。**

## ★機構 —— 2 つの相対和

底変換の係数は原文の慣習どおり

- 有限素点: `w | v` の係数は **`e(w|v)·c_v`**
- 無限素点: `w | v` の係数は **`[K_w : F_v]·c_v = (mult w / mult v)·c_v`**

で、次数が `[K:F]` 倍になることは、それぞれ

- `sum_ramification_inertia_hos` —— `Σ_{w|v} e(w|v)·f(w|v) = [K:F]`
  (`FinitePlaceRel.lean`。★mathlib の基本等式を `HeightOneSpectrum` に翻訳したもの)
- `sum_mult_comap_eq` —— `Σ_{w|v} mult w = mult v·[K:F]`
  (`InfinitePlaceRel.lean`。★mathlib は絶対版しか持たない)

に**帰着する**。

★有限側で `f(w|v)` が要るのは、`deg` の重みが `log q_w` であって
`log q_v` ではないからである——**`q_w = q_v^{f(w|v)}`**
(`Ideal.absNorm_eq_pow_inertiaDeg_of_liesOver`、mathlib にある)。

★★**`e` と `f` が両方要る**のはここに理由がある。係数から `e`、重みから `f` が出て、
基本等式の `e·f` がちょうど揃う。
-/

namespace ABC3.Found.GenEll

open NumberField IsDedekindDomain

section BaseChange

variable (F K : Type*) [Field F] [NumberField F] [Field K] [NumberField K] [Algebra F K]

/-! ## ★分岐指数と剰余次数 -/

/-- `w | v` の**分岐指数** `e(w|v)`。 -/
noncomputable def ramIdxOver (w : FinitePlace K) : ℕ :=
  (hosComap F K w).asIdeal.ramificationIdx w.asIdeal

/-- `w | v` の**剰余次数** `f(w|v)`。 -/
noncomputable def inertiaOver (w : FinitePlace K) : ℕ :=
  (hosComap F K w).asIdeal.inertiaDeg w.asIdeal

/-- ★**`q_w = q_v^{f(w|v)}`** の対数版。

★`Ideal.absNorm_eq_pow_inertiaDeg_of_liesOver` が mathlib に**ある**ので、
`LiesOver` インスタンス(`hosComap_liesOver`)を渡すだけで済む。 -/
theorem log_residueCard_over (w : FinitePlace K) :
    Real.log (residueCard w)
      = (inertiaOver F K w : ℝ) * Real.log (residueCard (hosComap F K w)) := by
  have h : Ideal.absNorm w.asIdeal
      = Ideal.absNorm (hosComap F K w).asIdeal ^ inertiaOver F K w :=
    Ideal.absNorm_eq_pow_inertiaDeg_of_liesOver w.asIdeal (hosComap F K w).asIdeal
      (hosComap F K w).isPrime (hosComap F K w).ne_bot
  simp only [residueCard, h]
  push_cast
  rw [Real.log_pow]

/-! ## ★底変換した算術因子 -/

open scoped Classical in
/-- ★底変換の**有限側** —— `w` の係数は `e(w|v)·c_v`(`v = hosComap w`)。

★台が有限であることは、`a` の台の各 `v` の上に有限個しか `w` が無いこと
(`hosFiber`)から出る。 -/
noncomputable def baseChangeFin (a : ADiv F) : FinitePlace K →₀ ℤ :=
  Finsupp.onFinset (a.fin.support.biUnion (fun v => hosFiber F K v))
    (fun w => (ramIdxOver F K w : ℤ) * a.fin (hosComap F K w))
    (by
      intro w hw
      rw [Finset.mem_biUnion]
      refine ⟨hosComap F K w, ?_, (mem_hosFiber F K).2 rfl⟩
      rw [Finsupp.mem_support_iff]
      intro hcon
      exact hw (by rw [hcon, mul_zero]))

open scoped Classical in
/-- ★底変換の**無限側** —— `w` の係数は `[K_w : F_v]·c_v = (mult w / mult v)·c_v`。

★`InfinitePlace K` は有限型なので、台の有限性は自明である。 -/
noncomputable def baseChangeArc (a : ADiv F) : InfinitePlace K →₀ ℝ :=
  Finsupp.onFinset Finset.univ
    (fun w => ((InfinitePlace.mult w : ℝ)
        / (InfinitePlace.mult (w.comap (algebraMap F K)) : ℝ))
      * a.arc (w.comap (algebraMap F K)))
    (fun w _ => Finset.mem_univ w)

open scoped Classical in
/-- ★★**算術因子の底変換** `ADiv(F) → ADiv(K)`。 -/
noncomputable def baseChange (a : ADiv F) : ADiv K :=
  (baseChangeFin F K a, baseChangeArc F K a)

open scoped Classical in
@[simp] theorem baseChange_fin (a : ADiv F) :
    (baseChange F K a).fin = baseChangeFin F K a := rfl

open scoped Classical in
@[simp] theorem baseChange_arc (a : ADiv F) :
    (baseChange F K a).arc = baseChangeArc F K a := rfl

/-! ## ★有限側の計算 -/

open scoped Classical in
/-- ★`v ≠ v'` なら `hosFiber v` と `hosFiber v'` は交わらない。 -/
theorem hosFiber_pairwiseDisjoint (s : Finset (FinitePlace F)) :
    (↑s : Set (FinitePlace F)).PairwiseDisjoint (fun v => hosFiber F K v) := by
  intro v _ v' _ hne
  simp only [Function.onFun, Finset.disjoint_left]
  intro w hw hw'
  rw [mem_hosFiber] at hw hw'
  exact hne (hw ▸ hw')

open scoped Classical in
/-- ★★**底変換で有限側の次数は `[K:F]` 倍になる**。

★係数から `e`、重み `log q_w` から `f` が出て、
基本等式 `Σ_{w|v} e·f = [K:F]` がちょうど揃う。 -/
theorem sum_fin_base_change (a : ADiv F) :
    (baseChangeFin F K a).sum (fun w n => (n : ℝ) * Real.log (residueCard w))
      = (Module.finrank F K : ℝ)
        * a.fin.sum (fun v n => (n : ℝ) * Real.log (residueCard v)) := by
  classical
  rw [baseChangeFin, Finsupp.onFinset_sum _ (fun w => by push_cast; ring)]
  rw [Finset.sum_biUnion (hosFiber_pairwiseDisjoint F K a.fin.support)]
  rw [Finsupp.sum, Finset.mul_sum]
  refine Finset.sum_congr rfl fun v _ => ?_
  -- fiber の各項で `hosComap w = v` を使う
  have hterm : ∀ w ∈ hosFiber F K v,
      (((ramIdxOver F K w : ℤ) * a.fin (hosComap F K w) : ℤ) : ℝ)
          * Real.log (residueCard w)
        = ((ramIdxOver F K w * inertiaOver F K w : ℕ) : ℝ)
            * (((a.fin v : ℤ) : ℝ) * Real.log (residueCard v)) := by
    intro w hw
    rw [mem_hosFiber] at hw
    rw [log_residueCard_over F K w, hw]
    push_cast
    ring
  rw [Finset.sum_congr rfl hterm, ← Finset.sum_mul]
  have hsum : ∑ w ∈ hosFiber F K v,
      ((ramIdxOver F K w * inertiaOver F K w : ℕ) : ℝ)
        = ((Module.finrank F K : ℕ) : ℝ) := by
    rw [← Nat.cast_sum]
    congr 1
    rw [← sum_ramification_inertia_hos F K v]
    refine Finset.sum_congr rfl fun w hw => ?_
    rw [mem_hosFiber] at hw
    rw [ramIdxOver, inertiaOver, hw]
  rw [hsum]

/-! ## ★無限側の計算 -/

open scoped Classical in
/-- ★★**底変換で無限側の次数も `[K:F]` 倍になる**。

★`InfinitePlaceRel.lean` の `sum_arc_base_change` そのものである。 -/
theorem sum_arc_base_change_adiv (a : ADiv F) :
    (baseChangeArc F K a).sum (fun _ r => r)
      = (Module.finrank F K : ℝ) * a.arc.sum (fun _ r => r) := by
  classical
  rw [baseChangeArc, Finsupp.onFinset_sum _ (fun _ => rfl)]
  rw [sum_arc_base_change F K (fun v => a.arc v)]
  congr 1
  rw [Finsupp.sum_fintype]
  intro _
  rfl

/-! ## ★★底変換不変性 -/

open scoped Classical in
/-- ★★**`deg_K(a_K) = [K:F]·deg_F(a)`**。

原文 (GenEll p.4):
> — where xF : Spec(OF ) → X is any morphism that gives rise to x.

★有限側(`sum_fin_base_change`)と無限側(`sum_arc_base_change_adiv`)を足すだけである。 -/
theorem deg_baseChange (a : ADiv F) :
    deg (baseChange F K a) = (Module.finrank F K : ℝ) * deg a := by
  classical
  simp only [deg, baseChange_fin, baseChange_arc, sum_fin_base_change,
    sum_arc_base_change_adiv]
  ring

open scoped Classical in
/-- ★★★**正規化次数は底変換で不変である**。

`degNormalized(a) ≝ deg_F(a)/[F:ℚ]` について

> **`degNormalized_K(a_K) = degNormalized_F(a)`**

★`ArithDiv.lean` の `degNormalized` の docstring は
「正規化する理由は**有限次拡大で不変にする**ためであり、
★**不変性そのものは本ファイルでは示していない**」と書いていた。
**本定理がそれを示す。**

★機構は `[K:ℚ] = [F:ℚ]·[K:F]`(`Module.finrank_mul_finrank`)で
`[K:F]` が約分されることだけである。 -/
theorem degNormalized_baseChange (a : ADiv F) :
    degNormalized (baseChange F K a) = degNormalized a := by
  classical
  haveI : FiniteDimensional ℚ F := inferInstance
  haveI : FiniteDimensional F K := Module.Finite.of_restrictScalars_finite ℚ F K
  have htower : Module.finrank ℚ F * Module.finrank F K = Module.finrank ℚ K :=
    Module.finrank_mul_finrank ℚ F K
  have hFK : (0 : ℝ) < (Module.finrank F K : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := F) (M := K)
  have hF : (0 : ℝ) < (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := ℚ) (M := F)
  rw [degNormalized, degNormalized, deg_baseChange, ← htower]
  push_cast
  field_simp

end BaseChange

/-! ## ★出典の紐付け(`.src`) -/

def deg_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(次数写像 deg_F)",
    sectionId := "genell-deg" }

def degNormalized_baseChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(正規化次数 deg_F)",
    sectionId := "genell-deg" }

end ABC3.Found.GenEll
