/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.HalfShiftProd
import ABC3.Meta.Claim

/-!
# 第 1399 ブロック —— **同種のノルムの積は代表系の置換で不変**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——**同種のノルム**の下ごしらえ

残る 1 本は**同種のノルム**

    F_i(z) ≔ ∏_{w ∈ T} (℘_Λ(z+w) − e_i) = c_i · (℘_{Λ′}(z) − e′_i),
    c_i = ∏_{w ∈ T∖{0}} (℘_Λ(w) − e_i)

である。☆Liouville（`elliptic_liouville`、第 598）で出すには

| 段 | 内容 | 状態 |
|---|---|---|
| (a) | `F_i` は `Λ′`-周期的 | ★**本ブロック** |
| (b) | `R ≔ ∏_{w∈T∖0}(℘(z+w) − e_i)` は**偶関数**（したがって `R′(0) = 0`） | ★**本ブロック** |
| (c) | `F_i − c_i·℘_{Λ′}` が原点で解析的（`R − c_i` が 2 位で消えることによる） | ☆次 |
| (d) | Liouville ＋ `F_i(v_i) = 0` | ☆その次 |

★★★(a)(b) はどちらも「代表系の置換」で出る——和の版（`veluAnalyticX_shift`・
`veluAnalyticX_periodic`、第 599・602）の**積版**である。

☆(b) の**負の代表**（`w ↦ −w` が誘導する `T` の置換）は新規に作った
——和の版では要らなかったものである。
-/

namespace ABC3.Found.GenEll

open PeriodPair Finset ABC3.Meta

open scoped Classical

/-! ## ★★★★負の代表 -/

/-- ★★★★**負の代表**——`w ↦ (−w の `T` における代表)`。 -/
noncomputable def negRep (P : PeriodPair) (T : Finset ℂ) (w : ℂ) : ℂ :=
  if h : ∃ v ∈ T, w + v ∈ P.lattice then h.choose else 0

/-- ★★★★★★★★★★★★
**負号は代表系を置換する**——★**代表系であることだけから**（第 1399）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆単射性は代表系の一意性から、全射性は `T` が有限だから出る
——`exists_veluShiftPerm`（第 602）と同じ型である。 -/
theorem exists_veluNegPerm (P P' : PeriodPair) (T : Finset ℂ)
    (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice) :
    ∃ σ : ℂ → ℂ, (∀ w ∈ T, σ w ∈ T) ∧ (∀ w ∈ T, ∀ w' ∈ T, σ w = σ w' → w = w')
      ∧ (∀ v ∈ T, ∃ w, ∃ _ : w ∈ T, σ w = v) ∧ (∀ w ∈ T, w + σ w ∈ P.lattice) := by
  classical
  have hex : ∀ w ∈ T, ∃ v ∈ T, w + v ∈ P.lattice := by
    intro w hw
    obtain ⟨v, hv, hv2, -⟩ := hrep w (hT w hw)
    exact ⟨v, hv, hv2⟩
  have hmemT : ∀ w ∈ T, negRep P T w ∈ T := by
    intro w hw
    rw [negRep, dif_pos (hex w hw)]
    exact (hex w hw).choose_spec.1
  have hprop : ∀ w ∈ T, w + negRep P T w ∈ P.lattice := by
    intro w hw
    rw [negRep, dif_pos (hex w hw)]
    exact (hex w hw).choose_spec.2
  have hinj : ∀ w ∈ T, ∀ w' ∈ T, negRep P T w = negRep P T w' → w = w' := by
    intro w hw w' hw' h
    have s1 := hprop w hw
    have s2 := hprop w' hw'
    rw [h] at s1
    have hdiff : -w' + w ∈ P.lattice := by
      have hs := P.lattice.sub_mem s1 s2
      have he : (w + negRep P T w') - (w' + negRep P T w') = -w' + w := by ring
      rw [he] at hs
      exact hs
    have hzero : -w' + w' ∈ P.lattice := by simp
    exact (rep_unique P P' T hrep (-w') (neg_mem (hT w' hw')) w' w hw' hw hzero hdiff).symm
  refine ⟨negRep P T, hmemT, hinj, ?_, hprop⟩
  intro v hv
  obtain ⟨a, ha, hav⟩ := Finset.surj_on_of_inj_on_of_card_le (s := T) (t := T)
    (fun a _ => negRep P T a) (fun a ha => hmemT a ha)
    (fun a₁ a₂ ha₁ ha₂ h => hinj a₁ ha₁ a₂ ha₂ h) le_rfl v hv
  exact ⟨a, ha, hav.symm⟩

/-! ## ★★★★同種のノルムの積 -/

/-- ★★★★★★**同種のノルムの積** `F(z) = ∏_{w∈T}(℘(z+w) − e)`。 -/
noncomputable def veluProd (P : PeriodPair) (T : Finset ℂ) (e z : ℂ) : ℂ :=
  ∏ w ∈ T, (P.weierstrassP (z + w) - e)

/-- ★★★★★★★★**平行移動で積は不変**（第 1399）。 -/
theorem veluProd_perm (P : PeriodPair) (T : Finset ℂ) (e a : ℂ) (σ : ℂ → ℂ)
    (hmem : ∀ w ∈ T, σ w ∈ T) (hinj : ∀ w ∈ T, ∀ w' ∈ T, σ w = σ w' → w = w')
    (hsurj : ∀ v ∈ T, ∃ w, ∃ _ : w ∈ T, σ w = v)
    (hshift : ∀ w ∈ T, a + w - σ w ∈ P.lattice) (z : ℂ) :
    ∏ w ∈ T, (P.weierstrassP (z + a + w) - e) = veluProd P T e z := by
  refine Finset.prod_nbij (i := σ) (fun w hw => hmem w hw)
    (fun w hw w' hw' h => hinj w hw w' hw' h) ?_ ?_
  · intro v hv
    obtain ⟨w, hw, hwv⟩ := hsurj v hv
    exact ⟨w, hw, hwv⟩
  · intro w hw
    have hd : z + a + w = (z + σ w) + (a + w - σ w) := by ring
    rw [hd, P.weierstrassP_add_coe _ ⟨_, hshift w hw⟩]

/-- ★★★★★★★★★★★★
**`F` は `Λ′`-周期的**——★**代表系であることだけから**（第 1399）。 -/
theorem veluProd_periodic (P P' : PeriodPair) (T : Finset ℂ)
    (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    (e z : ℂ) {l : ℂ} (hl : l ∈ P'.lattice) :
    veluProd P T e (z + l) = veluProd P T e z := by
  obtain ⟨σ, hmem, hinj, hsurj, hshift⟩ := exists_veluShiftPerm P P' T hT hrep l hl
  have hEq : veluProd P T e (z + l) = ∏ w ∈ T, (P.weierstrassP (z + l + w) - e) := by
    simp only [veluProd]
  rw [hEq]
  exact veluProd_perm P T e l σ hmem hinj hsurj (by
    intro w hw
    have h := hshift w hw
    have he : l + w - σ w = w + l - σ w := by ring
    rw [he]; exact h) z

/-- ★★★★★★**負の置換で積は不変**（第 1399）。 -/
theorem veluProd_negPerm (P : PeriodPair) (S : Finset ℂ) (e : ℂ) (σ : ℂ → ℂ)
    (hmem : ∀ w ∈ S, σ w ∈ S) (hinj : ∀ w ∈ S, ∀ w' ∈ S, σ w = σ w' → w = w')
    (hsurj : ∀ v ∈ S, ∃ w, ∃ _ : w ∈ S, σ w = v)
    (hneg : ∀ w ∈ S, w + σ w ∈ P.lattice) (z : ℂ) :
    ∏ w ∈ S, (P.weierstrassP (-z + w) - e) = ∏ w ∈ S, (P.weierstrassP (z + w) - e) := by
  refine Finset.prod_nbij (i := σ) (fun w hw => hmem w hw)
    (fun w hw w' hw' h => hinj w hw w' hw' h) ?_ ?_
  · intro v hv
    obtain ⟨w, hw, hwv⟩ := hsurj v hv
    exact ⟨w, hw, hwv⟩
  · intro w hw
    have h1 : P.weierstrassP (-z + w) = P.weierstrassP (z - w) := by
      rw [show z - w = -(-z + w) by ring, P.weierstrassP_neg]
    have hd : z - w = (z + σ w) + (-(w + σ w)) := by ring
    rw [h1, hd, P.weierstrassP_add_coe (z + σ w) ⟨-(w + σ w), neg_mem (hneg w hw)⟩]

/-- ☆負の代表は `0` を動かさない（第 1399）。 -/
theorem negPerm_zero (P P' : PeriodPair) (T : Finset ℂ) (h0T : (0:ℂ) ∈ T)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    {σ : ℂ → ℂ} (hmem : ∀ w ∈ T, σ w ∈ T) (hneg : ∀ w ∈ T, w + σ w ∈ P.lattice) :
    σ 0 = 0 := by
  have h1 : (0:ℂ) + σ 0 ∈ P.lattice := hneg 0 h0T
  have h2 : (0:ℂ) + 0 ∈ P.lattice := by simp
  exact rep_unique P P' T hrep 0 P'.lattice.zero_mem (σ 0) 0 (hmem 0 h0T) h0T h1 h2

/-- ★★★★★★★★★★★★★★★★
**`R(z) = ∏_{w∈T∖0}(℘(z+w) − e)` は偶関数**——★（第 1399）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これが `R′(0) = 0`、したがって `R − c` が原点で **2 位**で消えることの根拠である。
☆和の版（第 599・602）では要らなかった段である。 -/
theorem veluProd_erase_even (P P' : PeriodPair) (T : Finset ℂ) (h0T : (0:ℂ) ∈ T)
    (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice) (e z : ℂ) :
    ∏ w ∈ T.erase 0, (P.weierstrassP (-z + w) - e)
      = ∏ w ∈ T.erase 0, (P.weierstrassP (z + w) - e) := by
  obtain ⟨σ, hmem, hinj, hsurj, hneg⟩ := exists_veluNegPerm P P' T hT hrep
  have hσ0 : σ 0 = 0 := negPerm_zero P P' T h0T hrep hmem hneg
  refine veluProd_negPerm P (T.erase 0) e σ ?_ ?_ ?_ ?_ z
  · intro w hw
    refine Finset.mem_erase.2 ⟨?_, hmem w (Finset.mem_of_mem_erase hw)⟩
    intro hc
    exact (Finset.ne_of_mem_erase hw)
      (hinj w (Finset.mem_of_mem_erase hw) 0 h0T (by rw [hc, hσ0]))
  · intro w hw w' hw' h
    exact hinj w (Finset.mem_of_mem_erase hw) w' (Finset.mem_of_mem_erase hw') h
  · intro v hv
    obtain ⟨w, hw, hwv⟩ := hsurj v (Finset.mem_of_mem_erase hv)
    refine ⟨w, Finset.mem_erase.2 ⟨?_, hw⟩, hwv⟩
    intro hc
    exact (Finset.ne_of_mem_erase hv) (by rw [← hwv, hc, hσ0])
  · intro w hw
    exact hneg w (Finset.mem_of_mem_erase hw)

/-! ## ★出典の紐付け(`.src`) -/

def negRep.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(負の代表——−w の T における代表)",
    sectionId := "genell-lemma-3-5" }

def exists_veluNegPerm.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(負号は代表系を置換する。★代表系であることだけから)",
    sectionId := "genell-lemma-3-5" }

def veluProd.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(同種のノルムの積 F(z) = ∏(℘(z+w) − e))",
    sectionId := "genell-lemma-3-5" }

def veluProd_perm.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(平行移動で積は不変。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluProd_periodic.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(同種のノルムの積は Λ′-周期的。★代表系であることだけから)",
    sectionId := "genell-lemma-3-5" }

def veluProd_negPerm.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(負の置換で積は不変。★無条件)",
    sectionId := "genell-lemma-3-5" }

def negPerm_zero.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(負の代表は 0 を動かさない。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluProd_erase_even.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(R(z) = ∏_{w∈T∖0}(℘(z+w) − e) は偶関数。★代表系であることだけから)",
    sectionId := "genell-lemma-3-5" }

def veluProd_erase_even.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_veluShiftPerm(第 602、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_veluShiftPerm") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1399）**——同種のノルム " ++
       "`∏_{w∈T}(℘(z+w) − e_i) = c_i(℘_{Λ′}(z) − e′_i)` を Liouville で出すための下ごしらえ。" ++
       "☆`R` が偶であることから `R′(0) = 0`、したがって `R − c_i` は原点で 2 位で消え、" ++
       "`F_i − c_i℘_{Λ′}` の極が打ち消し合う。" ++
       "★和の版（第 599・602）では要らなかった段である。") 17 ]

end ABC3.Found.GenEll
