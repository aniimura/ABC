/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluProdPerm
import ABC3.Meta.Claim

/-!
# 第 1400 ブロック —— **同種のノルム**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★★★★★これは何か——**恒等式の心臓**

    ∏_{w ∈ T} (℘_Λ(z+w) − ℘_Λ(v)) = ( ∏_{w ∈ T∖{0}} (℘_Λ(w) − ℘_Λ(v)) ) · (℘_{Λ′}(z) − ℘_{Λ′}(v))

（`v ∉ Λ′`・`z ∉ Λ′`）。★これが `Δ(E)^l = Δ(E/C)·N⁴` の証明の**第 3 の部品**であり、
第 1397（判別式の差積表示）・第 1398（`∏_{i≠j}(℘(v_j+w)−e_i) = −D`）と合わせて
恒等式が出る。

☆古典的には「`x − e_i` の同種に沿ったノルムは `X − e′_i` に等しい」という主張であり、
`σ` 函数でも因子の理論でも書けるが、★**ここでは Liouville だけで出した**。

## ★★★★★★★★証明の骨

両辺とも `Λ′`-楕円で、`Λ′` の各点に 2 位の極、`v + Λ′` に 2 位の零点をもつ。
差 `W ≔ 左辺 − c·(℘_{Λ′} − e′)` は極が打ち消し合って整になり、`W(v) = 0` から `W ≡ 0`。

★★極の打ち消しの核心は **`R(z) ≔ ∏_{w∈T∖0}(℘(z+w) − e)` が原点で 2 位で `c` に近づく**
ことである:

* `R(0) = c`（定義そのもの）
* **`R` は偶関数**（第 1399——負の代表の置換）だから `R′(0) = 0`

☆したがって `R − c` の解析的位数は `≥ 2`（`natCast_le_analyticOrderAt_iff_iteratedDeriv_eq_zero`）、
すなわち `R − c = z²·g`（`g` 整）であり、`℘(z)·(R(z) − c) = g(z) + (℘(z) − 1/z²)(R(z) − c)` は整。

★★★**Lean 特有の注意**——`℘` は格子点で 0（ジャンク値）なので、
`W` そのものは格子点で連続にならない。☆そこで `Λ′` の上で値を `E(0)`（解析接続の値）に
置き換えた `f` に Liouville を当てる。`f` の格子点での解析性は**周期性で原点に帰着**する。
-/

namespace ABC3.Found.GenEll

open PeriodPair Finset ABC3.Meta Filter Topology

open scoped Classical

/-! ## ★★★★下ごしらえ -/

/-- ☆`0 ∈ T` なら積は `w = 0` の因子で割れる。 -/
theorem veluProd_split (P : PeriodPair) (T : Finset ℂ) (h0T : (0:ℂ) ∈ T) (e z : ℂ) :
    veluProd P T e z
      = (P.weierstrassP z - e) * ∏ w ∈ T.erase 0, (P.weierstrassP (z + w) - e) := by
  rw [veluProd, ← Finset.mul_prod_erase _ _ h0T, add_zero]

/-- ☆`℘_Λ − ℘_{Λ′}` は原点で解析的——主要部が同じだから。 -/
theorem weierstrassP_sub_analyticAt_zero (P P' : PeriodPair) :
    AnalyticAt ℂ (fun z => P.weierstrassP z - P'.weierstrassP z) 0 := by
  have h1 : AnalyticAt ℂ (fun z => P.weierstrassP z - 1 / (z - 0) ^ 2) 0 :=
    weierstrassP_sub_pole_analyticAt P 0 P.lattice.zero_mem
  have h2 : AnalyticAt ℂ (fun z => P'.weierstrassP z - 1 / (z - 0) ^ 2) 0 :=
    weierstrassP_sub_pole_analyticAt P' 0 P'.lattice.zero_mem
  refine (h1.sub h2).congr ?_
  filter_upwards with z
  simp only [Pi.sub_apply]
  ring

/-- ☆`Λ′` の外では積は解析的。 -/
theorem veluProd_analyticAt_of_notMem (P P' : PeriodPair) (hle : P.lattice ≤ P'.lattice)
    (T : Finset ℂ) (hT : ∀ w ∈ T, w ∈ P'.lattice) (e : ℂ) {p : ℂ} (hp : p ∉ P'.lattice) :
    AnalyticAt ℂ (veluProd P T e) p := by
  refine Finset.analyticAt_fun_prod _ fun w hw => ?_
  refine (shifted_analyticAt P p w ?_).sub analyticAt_const
  exact notMem_of_notMem P P' hle T hT p hp w hw

/-- ☆偶関数の原点での微分は `0`。 -/
theorem deriv_eq_zero_of_even (f : ℂ → ℂ) (h : ∀ z, f (-z) = f z) : deriv f 0 = 0 := by
  have h1 : deriv (fun z => f (-z)) 0 = -deriv f 0 := by
    simpa using (deriv_comp_neg f 0)
  have h2 : (fun z => f (-z)) = f := funext h
  rw [h2] at h1
  linear_combination h1 / 2

/-- ★★★★★★**偶関数は定数値に 2 位で近づく**（第 1400）。 -/
theorem two_le_order_of_even {f : ℂ → ℂ} {c : ℂ} (hana : AnalyticAt ℂ (fun z => f z - c) 0)
    (h0 : f 0 = c) (heven : ∀ z, f (-z) = f z) :
    ((2 : ℕ) : ℕ∞) ≤ analyticOrderAt (fun z => f z - c) 0 := by
  rw [natCast_le_analyticOrderAt_iff_iteratedDeriv_eq_zero hana]
  intro i hi
  interval_cases i
  · simp only [iteratedDeriv_zero]; rw [h0]; ring
  · rw [iteratedDeriv_one]
    exact deriv_eq_zero_of_even _ (fun z => by rw [heven z])

/-- ☆`R` は原点で解析的。 -/
theorem veluProdErase_analyticAt (P P' : PeriodPair) (T : Finset ℂ) (h0T : (0:ℂ) ∈ T)
    (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice) (e : ℂ) :
    AnalyticAt ℂ (fun z => ∏ w ∈ T.erase 0, (P.weierstrassP (z + w) - e)) 0 := by
  refine Finset.analyticAt_fun_prod _ fun w hw => ?_
  refine (shifted_analyticAt P 0 w ?_).sub analyticAt_const
  simpa using rep_notMem_lattice P P' T h0T hT hrep hw

/-- ☆`R(0) = c`。 -/
theorem veluProdErase_zero (P : PeriodPair) (T : Finset ℂ) (e : ℂ) :
    (∏ w ∈ T.erase 0, (P.weierstrassP ((0:ℂ) + w) - e))
      = ∏ w ∈ T.erase 0, (P.weierstrassP w - e) := by
  refine Finset.prod_congr rfl fun w _ => ?_
  rw [zero_add]

/-- ★★★★★★★★**補正した函数は整**——格子点では周期性で原点に帰着する（第 1400）。 -/
theorem veluNorm_differentiable (P P' : PeriodPair) (hle : P.lattice ≤ P'.lattice)
    (T : Finset ℂ) (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (e e' c : ℂ) (E : ℂ → ℂ)
    (W : ℂ → ℂ) (hW : W = fun y => veluProd P T e y - c * (P'.weierstrassP y - e'))
    (f : ℂ → ℂ) (hf : f = fun y => if y ∈ P'.lattice then E 0 else W y)
    (hfper : ∀ y : ℂ, ∀ l ∈ P'.lattice, f (y + l) = f y)
    (hf0 : AnalyticAt ℂ f 0) :
    Differentiable ℂ f := by
  intro p
  by_cases hp : p ∈ P'.lattice
  · have hshift : (fun y => f (y - p)) = f := by
      funext y
      have hh := hfper (y - p) p hp
      rw [sub_add_cancel] at hh
      exact hh.symm
    have hcomp : AnalyticAt ℂ (fun y : ℂ => f (y - p)) p := by
      refine AnalyticAt.comp (g := f) (f := fun y : ℂ => y - p) (x := p) ?_
        (analyticAt_id.sub analyticAt_const)
      simpa using hf0
    rw [hshift] at hcomp
    exact hcomp.differentiableAt
  · have hopen : IsOpen ((P'.lattice : Set ℂ)ᶜ) := P'.isClosed_lattice.isOpen_compl
    have hWana : AnalyticAt ℂ W p := by
      rw [hW]
      exact (veluProd_analyticAt_of_notMem P P' hle T hT e hp).sub
        (analyticAt_const.mul ((P'.analyticOnNhd_weierstrassP p hp).sub analyticAt_const))
    refine (hWana.congr ?_).differentiableAt
    filter_upwards [hopen.mem_nhds hp] with y hy
    simp only [Set.mem_compl_iff, SetLike.mem_coe] at hy
    rw [hf]
    simp only []
    rw [if_neg hy]

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★同種のノルム -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] 同種のノルム**——★**代表系であることだけから**（第 1400）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

    ∏_{w∈T} (℘_Λ(z+w) − ℘_Λ(v)) = ( ∏_{w∈T∖0} (℘_Λ(w) − ℘_Λ(v)) ) · (℘_{Λ′}(z) − ℘_{Λ′}(v))

★★★これが `Δ(E)^l = Δ(E/C)·N⁴` の**心臓**である。
☆`σ` 函数も因子の理論も要らず、**Liouville（第 598）と偶関数性（第 1399）だけ**で出た。 -/
theorem veluProd_eq (P P' : PeriodPair) (hle : P.lattice ≤ P'.lattice)
    (T : Finset ℂ) (h0T : (0:ℂ) ∈ T) (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    {v : ℂ} (hv : v ∉ P'.lattice) {z : ℂ} (hz : z ∉ P'.lattice) :
    veluProd P T (P.weierstrassP v) z
      = (∏ w ∈ T.erase 0, (P.weierstrassP w - P.weierstrassP v))
        * (P'.weierstrassP z - P'.weierstrassP v) := by
  classical
  set e := P.weierstrassP v with he
  set c := ∏ w ∈ T.erase 0, (P.weierstrassP w - e) with hc
  set e' := P'.weierstrassP v with he'
  set R : ℂ → ℂ := fun y => ∏ w ∈ T.erase 0, (P.weierstrassP (y + w) - e) with hR
  set W : ℂ → ℂ := fun y => veluProd P T e y - c * (P'.weierstrassP y - e') with hW
  have hRana : AnalyticAt ℂ R 0 := veluProdErase_analyticAt P P' T h0T hT hrep e
  have hR0 : R 0 = c := veluProdErase_zero P T e
  have hReven : ∀ y, R (-y) = R y := fun y => veluProd_erase_even P P' T h0T hT hrep e y
  have h2le := two_le_order_of_even (hRana.sub analyticAt_const) hR0 hReven
  obtain ⟨g, hgana, hgeq⟩ := (natCast_le_analyticOrderAt (hRana.sub analyticAt_const)).mp h2le
  have hgeq' : ∀ᶠ y in 𝓝 (0:ℂ), R y - c = y ^ 2 * g y := by
    filter_upwards [hgeq] with y hy
    simpa [smul_eq_mul] using hy
  set A : ℂ → ℂ := fun y => P.weierstrassP y - 1 / y ^ 2 with hA
  set E : ℂ → ℂ := fun y => c * (P.weierstrassP y - P'.weierstrassP y) - e * R y + c * e'
      + g y + A y * (R y - c) with hE
  have hWE : ∀ y : ℂ, y ≠ 0 → R y - c = y ^ 2 * g y → W y = E y := by
    intro y hy0 hfac
    have hsplit : veluProd P T e y = (P.weierstrassP y - e) * R y :=
      veluProd_split P T h0T e y
    have h1 : R y = c + y ^ 2 * g y := by linear_combination hfac
    simp only [hW, hE, hA]
    rw [hsplit, h1]
    field_simp
    ring
  have hWper : ∀ y : ℂ, ∀ l ∈ P'.lattice, W (y + l) = W y := by
    intro y l hl
    simp only [hW]
    rw [veluProd_periodic P P' T hT hrep e y hl, P'.weierstrassP_add_coe y ⟨l, hl⟩]
  set f : ℂ → ℂ := fun y => if y ∈ P'.lattice then E 0 else W y with hf
  have hfper : ∀ y : ℂ, ∀ l ∈ P'.lattice, f (y + l) = f y := by
    intro y l hl
    simp only [hf]
    by_cases hy : y ∈ P'.lattice
    · rw [if_pos (P'.lattice.add_mem hy hl), if_pos hy]
    · have hyl : y + l ∉ P'.lattice := fun hcc => hy (by simpa using P'.lattice.sub_mem hcc hl)
      rw [if_neg hyl, if_neg hy]
      exact hWper y l hl
  have hAana : AnalyticAt ℂ A 0 := by
    refine (weierstrassP_sub_pole_analyticAt P 0 P.lattice.zero_mem).congr ?_
    filter_upwards with y
    simp only [hA, sub_zero]
  have hEana : AnalyticAt ℂ E 0 := by
    simp only [hE]
    exact ((((analyticAt_const.mul (weierstrassP_sub_analyticAt_zero P P')).sub
      (analyticAt_const.mul hRana)).add analyticAt_const).add hgana).add
      (hAana.mul (hRana.sub analyticAt_const))
  have hnbhd : ((P'.lattice : Set ℂ) \ {0})ᶜ ∈ 𝓝 (0:ℂ) :=
    P'.isOpen_compl_lattice_sdiff.mem_nhds (by simp)
  have hfE : f =ᶠ[𝓝 (0:ℂ)] E := by
    filter_upwards [hnbhd, hgeq'] with y hy hfac
    by_cases hy0 : y = 0
    · subst hy0; simp only [hf]; rw [if_pos P'.lattice.zero_mem]
    · have hyΛ' : y ∉ P'.lattice := fun hcc => hy ⟨hcc, hy0⟩
      simp only [hf]
      rw [if_neg hyΛ']
      exact hWE y hy0 hfac
  have hf0 : AnalyticAt ℂ f 0 := hEana.congr hfE.symm
  have hfdiff : Differentiable ℂ f :=
    veluNorm_differentiable P P' hle T hT e e' c E W hW f hf hfper hf0
  have hall := elliptic_liouville P' f hfdiff hfper
  have hveluv : veluProd P T e v = 0 := by
    rw [veluProd_split P T h0T e v, he]
    simp
  have hfv : f v = 0 := by
    simp only [hf]
    rw [if_neg hv]
    simp only [hW]
    rw [hveluv, he']
    ring
  have hfz := hall z v
  rw [hfv] at hfz
  simp only [hf] at hfz
  rw [if_neg hz] at hfz
  simp only [hW] at hfz
  linear_combination hfz

/-! ## ★出典の紐付け(`.src`) -/

def veluProd_split.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(0 ∈ T なら積は w = 0 の因子で割れる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def weierstrassP_sub_analyticAt_zero.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(℘_Λ − ℘_{Λ′} は原点で解析的。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluProd_analyticAt_of_notMem.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Λ′ の外では積は解析的。★無条件)",
    sectionId := "genell-lemma-3-5" }

def deriv_eq_zero_of_even.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(偶関数の原点での微分は 0。★無条件)",
    sectionId := "genell-lemma-3-5" }

def two_le_order_of_even.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(偶関数は定数値に 2 位で近づく。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluProdErase_analyticAt.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(R は原点で解析的。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluProdErase_zero.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(R(0) = c。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluNorm_differentiable.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(補正した函数は整——格子点では周期性で原点に帰着。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluProd_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(同種のノルム。★代表系であることだけから)",
    sectionId := "genell-lemma-3-5" }

def veluProd_eq.needs : List ProofObligation :=
  [ .citation "[ABC3]" "elliptic_liouville(第 598、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.elliptic_liouville") 1,
    .citation "[ABC3]" "veluProd_erase_even(第 1399、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.veluProd_erase_even") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1400）**——`Δ(E)^l = Δ(E/C)·N⁴` の心臓である。" ++
       "☆古典的には `σ` 函数か因子の理論で書く所を、" ++
       "**Liouville と `R` の偶関数性だけ**で出した。" ++
       "★極の打ち消しの核心は `R − c` の解析的位数が `≥ 2` であること" ++
       "（`R(0) = c` ＋ `R` が偶 ⟹ `R′(0) = 0`）。" ++
       "★★`℘` が格子点でジャンク値 `0` を取るので `W` は格子点で連続にならない" ++
       "——`Λ′` の上で解析接続の値に置き換えた `f` に Liouville を当てた。") 17 ]

end ABC3.Found.GenEll
