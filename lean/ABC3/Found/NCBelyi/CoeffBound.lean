import ABC3.Meta.Claim
import Mathlib.RingTheory.Polynomial.Vieta
import Mathlib.Analysis.Complex.Polynomial.Basic

/-!
# [NCBelyi] Lemma 2.4 —— 最小多項式の係数と値の評価(`Found`)

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.5。

原文 (NCBelyi p.5):
> Then one verifies immediately that all of the coefficients of f0(x) have absolute value

## ★★本ファイルが取るもの

`Lemma 2.4` の帰納段は、`α₀` の**モニックな最小多項式 `f₀`** を作って
`S′ ≝ f₀(S) ∪ f₀(S₀)` へ移る(`S₀` は `f₀′` の根の集合)。
そこで原文が畳んでいる定量的な段は 3 つある:

| 原文 | ここでの形 |
|---|---|
| 係数の絶対値は `≤ d₀^{d₀}` | `norm_coeff_le_choose` —— `≤ C(d₀,k)` |
| `f₀` の値は「`d₀` のある式で抑えられる」 | `norm_eval_le_pow` —— `≤ (‖x‖+1)^{d₀}` |
| 「`C` を適当に大きく取れば」`f₀(β) ∉ S′` | `norm_eval_lt_norm_eval_of_le` |

## ★★★逸脱の記録(強化、2026-08-27)

原文は係数の限界を **`d₀^{d₀}`** と書く。★我々が取るのは **`C(d₀,k) ≤ 2^{d₀}`** で、
`d₀ ≥ 2` では原文より**強い**。★★後続(値の評価と `C` の選択)は上界しか使わないので、
強い方を取っても影響しない。★★★原文の `d₀^{d₀}` は「粗くてよい」ことの表明であって、
主張の内容ではない——**粗さを再現する必要はない**。

## ★★★★根はすべて単位円板の中、という仮定について

原文は正規化(`Lemma 2.3` の自己同型 + 正の有理数倍)で
`|α| ≤ 1`(`∀ α ∈ S\{∞}`)を達成する。`S` は `Gal(ℚ̄/ℚ)`-安定に取ってあるので、
★**`α₀` の共役はすべて `S` に入る**——ゆえに `f₀` の根はすべて単位円板の中である。
本ファイルはその状況を `∀ z ∈ f.roots, ‖z‖ ≤ 1` として仮定に置く。

## ★★★★★上界は係数を経由しなくても出る

原文は「係数が `≤ d₀^{d₀}` ⇒ 値が抑えられる」と**係数を経由する**。
★実装では `f = ∏(X − αᵢ)` から直接 `‖f(x)‖ ≤ (‖x‖+1)^{d₀}` が出た(`norm_eval_le_pow`)
——係数の限界は要らない。★★それでも `norm_coeff_le_choose` を取るのは、
**`f₀′` の根の位置を Cauchy 限界で抑える**のに係数が要るからである(次のブロック)。
-/

namespace ABC3.Found.NCBelyi

open Polynomial

/-! ## ★単位円板の中の積 -/

/-- ★**単位円板の中の元の積は単位円板の中**。 -/
theorem norm_prod_le_one (t : Multiset ℂ) (h : ∀ z ∈ t, ‖z‖ ≤ 1) : ‖t.prod‖ ≤ 1 := by
  induction t using Multiset.induction with
  | empty => simp
  | cons a s ih =>
    rw [Multiset.prod_cons, norm_mul]
    have ha : ‖a‖ ≤ 1 := h a (Multiset.mem_cons_self a s)
    have hs : ‖s.prod‖ ≤ 1 := ih (fun z hz => h z (Multiset.mem_cons_of_mem hz))
    nlinarith [norm_nonneg a, norm_nonneg s.prod]

/-- ★★**基本対称式の限界** —— `|e_k| ≤ C(n,k)`。

★`esymm` は「大きさ `k` の部分多重集合の積の総和」なので、
項数が `C(n,k)` 個、各項が `≤ 1` である。それだけ。 -/
theorem norm_esymm_le (s : Multiset ℂ) (h : ∀ z ∈ s, ‖z‖ ≤ 1) (k : ℕ) :
    ‖s.esymm k‖ ≤ (s.card.choose k : ℝ) := by
  rw [Multiset.esymm]
  refine le_trans (norm_multiset_sum_le _) ?_
  rw [Multiset.map_map]
  have hle : ((Multiset.powersetCard k s).map (fun t => ‖t.prod‖)).sum
      ≤ ((Multiset.powersetCard k s).map (fun _ => (1 : ℝ))).sum := by
    refine Multiset.sum_map_le_sum_map _ _ ?_
    intro t ht
    exact norm_prod_le_one t
      (fun z hz => h z (Multiset.mem_of_le (Multiset.mem_powersetCard.1 ht).1 hz))
  refine le_trans hle ?_
  rw [Multiset.map_const', Multiset.sum_replicate, Multiset.card_powersetCard]
  simp

/-! ## ★★★★★原文の「係数の絶対値は `≤ d₀^{d₀}`」 -/

/-- ★★★★**係数の限界**。

原文 (NCBelyi p.5):
> Then one verifies immediately that all of the coefficients of f0(x) have absolute value

★原文の限界は `d₀^{d₀}`、ここでの限界は `C(d₀,k)`(★`≤ 2^{d₀}` なので原文より強い)。
★★Vieta(`Multiset.prod_X_sub_C_coeff`)で係数を基本対称式に書き換えるだけである。 -/
theorem norm_coeff_le_choose (f : ℂ[X]) (hf : f.Monic)
    (h : ∀ z ∈ f.roots, ‖z‖ ≤ 1) (k : ℕ) :
    ‖f.coeff k‖ ≤ (f.natDegree.choose k : ℝ) := by
  have hcard : f.roots.card = f.natDegree :=
    Polynomial.splits_iff_card_roots.1 (IsAlgClosed.splits f)
  by_cases hk : k ≤ f.natDegree
  · have hprod : (Multiset.map (fun a => X - C a) f.roots).prod = f :=
      Polynomial.prod_multiset_X_sub_C_of_monic_of_roots_card_eq hf hcard
    have hco := Multiset.prod_X_sub_C_coeff f.roots (k := k) (by omega)
    rw [hprod] at hco
    rw [hco, norm_mul, norm_pow, norm_neg, norm_one, one_pow, one_mul]
    refine le_trans (norm_esymm_le f.roots h _) ?_
    rw [hcard, Nat.choose_symm hk]
  · rw [Polynomial.coeff_eq_zero_of_natDegree_lt (by omega),
      Nat.choose_eq_zero_of_lt (by omega)]
    simp

/-! ## ★★★★★★値の上界と下界 -/

/-- ★**上界(多重集合版)** —— `‖∏(x − αᵢ)‖ ≤ (‖x‖+1)^n`。 -/
theorem norm_prod_eval_le (s : Multiset ℂ) (h : ∀ z ∈ s, ‖z‖ ≤ 1) (x : ℂ) :
    ‖(Multiset.map (fun a => X - C a) s).prod.eval x‖ ≤ (‖x‖ + 1) ^ s.card := by
  induction s using Multiset.induction with
  | empty => simp
  | cons a t ih =>
    have ha : ‖a‖ ≤ 1 := h a (Multiset.mem_cons_self a t)
    have iht := ih (fun z hz => h z (Multiset.mem_cons_of_mem hz))
    have hxa : ‖(X - C a : ℂ[X]).eval x‖ ≤ ‖x‖ + 1 := by
      rw [Polynomial.eval_sub, Polynomial.eval_X, Polynomial.eval_C]
      have := norm_sub_le x a
      linarith
    rw [Multiset.map_cons, Multiset.prod_cons, Polynomial.eval_mul, norm_mul,
      Multiset.card_cons, pow_succ']
    exact mul_le_mul hxa iht (norm_nonneg _) (by positivity)

/-- ★★★★**値の上界** —— 原文の「`d₀` のある式で抑えられる」の中身。

原文 (NCBelyi p.5):
> bounded by some fixed expression in d0. Thus, for a suitable choice of C, it follows

★★**係数を経由しない**——`f = ∏(X − αᵢ)` から直接出る。 -/
theorem norm_eval_le_pow (f : ℂ[X]) (hf : f.Monic) (h : ∀ z ∈ f.roots, ‖z‖ ≤ 1) (x : ℂ) :
    ‖f.eval x‖ ≤ (‖x‖ + 1) ^ f.natDegree := by
  have hcard : f.roots.card = f.natDegree :=
    Polynomial.splits_iff_card_roots.1 (IsAlgClosed.splits f)
  have hprod : (Multiset.map (fun a => X - C a) f.roots).prod = f :=
    Polynomial.prod_multiset_X_sub_C_of_monic_of_roots_card_eq hf hcard
  have hle := norm_prod_eval_le f.roots h x
  rw [hprod, hcard] at hle
  exact hle

/-- ★**単位円板の上では `≤ 2^{d₀}`** —— 上界の一番よく使う形。 -/
theorem norm_eval_le_two_pow (f : ℂ[X]) (hf : f.Monic)
    (h : ∀ z ∈ f.roots, ‖z‖ ≤ 1) {x : ℂ} (hx : ‖x‖ ≤ 1) :
    ‖f.eval x‖ ≤ 2 ^ f.natDegree := by
  refine le_trans (norm_eval_le_pow f hf h x) ?_
  refine pow_le_pow_left₀ (by positivity) (by linarith) _

/-- ★**下界(多重集合版)** —— `‖∏(x − αᵢ)‖ ≥ (‖x‖−1)^n`。 -/
theorem le_norm_prod_eval (s : Multiset ℂ) (h : ∀ z ∈ s, ‖z‖ ≤ 1) {x : ℂ} (hx : 1 ≤ ‖x‖) :
    (‖x‖ - 1) ^ s.card ≤ ‖(Multiset.map (fun a => X - C a) s).prod.eval x‖ := by
  induction s using Multiset.induction with
  | empty => simp
  | cons a t ih =>
    have ha : ‖a‖ ≤ 1 := h a (Multiset.mem_cons_self a t)
    have iht := ih (fun z hz => h z (Multiset.mem_cons_of_mem hz))
    have hxa : ‖x‖ - 1 ≤ ‖(X - C a : ℂ[X]).eval x‖ := by
      rw [Polynomial.eval_sub, Polynomial.eval_X, Polynomial.eval_C]
      have := norm_sub_norm_le x a
      linarith
    have h0 : (0 : ℝ) ≤ ‖x‖ - 1 := by linarith
    rw [Multiset.map_cons, Multiset.prod_cons, Polynomial.eval_mul, norm_mul,
      Multiset.card_cons, pow_succ']
    exact mul_le_mul hxa iht (pow_nonneg h0 _) (norm_nonneg _)

/-- ★★★**値の下界** —— `‖β‖` が大きければ `‖f₀(β)‖` も大きい。

★原文は書いていない段である。★★「`C` を適当に大きく取れば」を成り立たせているのは
**これ**である——上界だけでは `f₀(β) ∉ S′` は出ない。 -/
theorem le_norm_eval_of_one_le (f : ℂ[X]) (hf : f.Monic)
    (h : ∀ z ∈ f.roots, ‖z‖ ≤ 1) {x : ℂ} (hx : 1 ≤ ‖x‖) :
    (‖x‖ - 1) ^ f.natDegree ≤ ‖f.eval x‖ := by
  have hcard : f.roots.card = f.natDegree :=
    Polynomial.splits_iff_card_roots.1 (IsAlgClosed.splits f)
  have hprod : (Multiset.map (fun a => X - C a) f.roots).prod = f :=
    Polynomial.prod_multiset_X_sub_C_of_monic_of_roots_card_eq hf hcard
  have hle := le_norm_prod_eval f.roots h hx
  rw [hprod, hcard] at hle
  exact hle

/-! ## ★★★★★★★「suitable choice of C」 -/

/-- ★★★★★**`C` の選び方**。

原文 (NCBelyi p.5):
> bounded by some fixed expression in d0. Thus, for a suitable choice of C, it follows

`‖x‖ ≤ R` なるすべての `x` に対して `‖f₀(x)‖ < ‖f₀(β)‖` となるには、
★**`‖β‖ > R + 2` で足りる**。

★★原文は「`d₀` のある式」「適当な `C`」と 2 度畳むが、実測すると
`R + 2` という**`d₀` に依らない**形になった——`f₀` がモニックで根が単位円板の中だから、
`(R+1)^{d₀} < (‖β‖−1)^{d₀}` が `d₀` に関して一様に効く。
★★★これは原文より強い(原文は「`d₀` に相対的に十分大きい `C`」と書く)。

★★★★これが `f₀(β) ∉ f₀(S) ∪ f₀(S₀)` の中身である
——`S` は単位円板の中(`R = 1`)、`S₀` は `f₀′` の根(`R` は Cauchy 限界)。 -/
theorem norm_eval_lt_norm_eval_of_le (f : ℂ[X]) (hf : f.Monic) (hdeg : 0 < f.natDegree)
    (h : ∀ z ∈ f.roots, ‖z‖ ≤ 1) {x β : ℂ} (R : ℝ) (hx : ‖x‖ ≤ R) (hR : 0 ≤ R)
    (hβ : R + 2 < ‖β‖) :
    ‖f.eval x‖ < ‖f.eval β‖ := by
  have h1 : ‖f.eval x‖ ≤ (‖x‖ + 1) ^ f.natDegree := norm_eval_le_pow f hf h x
  have h1' : (‖x‖ + 1) ^ f.natDegree ≤ (R + 1) ^ f.natDegree :=
    pow_le_pow_left₀ (by positivity) (by linarith) _
  have h2 : (‖β‖ - 1) ^ f.natDegree ≤ ‖f.eval β‖ :=
    le_norm_eval_of_one_le f hf h (by linarith)
  have h3 : (R + 1) ^ f.natDegree < (‖β‖ - 1) ^ f.natDegree :=
    pow_lt_pow_left₀ (by linarith) (by linarith) (by omega)
  linarith

/-- ★★★★★★**分離の形** —— `f₀(β)` は `f₀` の「小さい点」の像には入らない。

★これが `Lemma 2.4` の証明の `f0(β) ∉ S′` そのものである。 -/
theorem eval_ne_of_norm_le (f : ℂ[X]) (hf : f.Monic) (hdeg : 0 < f.natDegree)
    (h : ∀ z ∈ f.roots, ‖z‖ ≤ 1) {x β : ℂ} (R : ℝ) (hx : ‖x‖ ≤ R) (hR : 0 ≤ R)
    (hβ : R + 2 < ‖β‖) :
    f.eval β ≠ f.eval x := by
  intro hcon
  have := norm_eval_lt_norm_eval_of_le f hf hdeg h R hx hR hβ
  rw [hcon] at this
  exact lt_irrefl _ this

/-! ## ★出典の紐付け(`.src`) -/

def norm_coeff_le_choose.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(f₀ の係数の絶対値の限界)",
    sectionId := "ncbelyi-lemma-2-4" }

def norm_eval_le_pow.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(f₀ の値が d₀ のある式で抑えられること)",
    sectionId := "ncbelyi-lemma-2-4" }

def norm_eval_lt_norm_eval_of_le.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(suitable choice of C)",
    sectionId := "ncbelyi-lemma-2-4" }

def eval_ne_of_norm_le.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(f₀(β) ∉ S′)",
    sectionId := "ncbelyi-lemma-2-4" }

end ABC3.Found.NCBelyi
