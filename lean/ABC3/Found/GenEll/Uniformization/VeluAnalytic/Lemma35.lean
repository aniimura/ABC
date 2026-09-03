/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Analysis.SpecialFunctions.Elliptic.Weierstrass
import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Basic
import ABC3.Found.GenEll.LatticeCurve
import ABC3.Found.GenEll.WeierstrassODE
import ABC3.Found.GenEll.Velu
import ABC3.Found.GenEll.PointVariableChange
import ABC3.Meta.Claim
import ABC3.Found.GenEll.Uniformization.Basic
import ABC3.Found.GenEll.Uniformization.VeluAnalytic.Proposition34

/-!
# VeluAnalytic —— `[GenEll] Lemma 3.5` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve PeriodPair

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★Vélu の公式の解析側 -/

/-- ★★★★★★**Vélu の `X` の解析側**

    `X(z) = Σ_{w ∈ T} ℘_Λ(z + w) − c`

★`T` は `Λ′/Λ` の代表系（`0` を含む）、`c = Σ_{w ∈ T∖{0}} ℘_Λ(w)` のつもりである。
原文の形 `℘(z) + Σ_{w≠0}[℘(z+w) − ℘(w)]` を、和をひとまとめにして書いたものである。 -/
noncomputable def veluAnalyticX (P : PeriodPair) (T : Finset ℂ) (c : ℂ) (z : ℂ) : ℂ :=
  (∑ w ∈ T, P.weierstrassP (z + w)) - c

/-- ★★★★★★★★★★★★**代表系は平行移動で置換される**——`℘` 側。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`σ` は「`w ↦ w + w₀` が `Λ′/Λ` に誘導する置換」を代表系の上に持ち上げたものである
（`hshift`: `w + w₀ − σ(w) ∈ Λ`）。 -/
theorem veluAnalyticSum_shift (P : PeriodPair) (T : Finset ℂ) (w₀ : ℂ) (σ : ℂ → ℂ)
    (hmem : ∀ w ∈ T, σ w ∈ T)
    (hinj : ∀ w ∈ T, ∀ w' ∈ T, σ w = σ w' → w = w')
    (hsurj : ∀ v ∈ T, ∃ w, ∃ _ : w ∈ T, σ w = v)
    (hshift : ∀ w ∈ T, w + w₀ - σ w ∈ P.lattice) (z : ℂ) :
    ∑ w ∈ T, P.weierstrassP (z + w + w₀) = ∑ w ∈ T, P.weierstrassP (z + w) := by
  refine Finset.sum_bij (fun w _ => σ w) (fun w hw => hmem w hw)
    (fun w hw w' hw' h => hinj w hw w' hw' h) hsurj ?_
  intro w hw
  have hz : z + w + w₀ = (z + σ w) + (w + w₀ - σ w) := by ring
  rw [hz]
  exact P.weierstrassP_add_coe _ ⟨_, hshift w hw⟩

/-- ★★★★★★★★★★★★**代表系は平行移動で置換される**——`℘′` 側。 -/
theorem derivVeluAnalyticSum_shift (P : PeriodPair) (T : Finset ℂ) (w₀ : ℂ) (σ : ℂ → ℂ)
    (hmem : ∀ w ∈ T, σ w ∈ T)
    (hinj : ∀ w ∈ T, ∀ w' ∈ T, σ w = σ w' → w = w')
    (hsurj : ∀ v ∈ T, ∃ w, ∃ _ : w ∈ T, σ w = v)
    (hshift : ∀ w ∈ T, w + w₀ - σ w ∈ P.lattice) (z : ℂ) :
    ∑ w ∈ T, P.derivWeierstrassP (z + w + w₀) = ∑ w ∈ T, P.derivWeierstrassP (z + w) := by
  refine Finset.sum_bij (fun w _ => σ w) (fun w hw => hmem w hw)
    (fun w hw w' hw' h => hinj w hw w' hw' h) hsurj ?_
  intro w hw
  have hz : z + w + w₀ = (z + σ w) + (w + w₀ - σ w) := by ring
  rw [hz]
  exact P.derivWeierstrassP_add_coe _ ⟨_, hshift w hw⟩

/-- ★★★★★★★★★★★★★★★★**Vélu の `X` は `Λ′`-周期的**。

★これが「`X` が `E/H` の座標である」ことの中身である。 -/
theorem veluAnalyticX_shift (P : PeriodPair) (T : Finset ℂ) (c w₀ : ℂ) (σ : ℂ → ℂ)
    (hmem : ∀ w ∈ T, σ w ∈ T)
    (hinj : ∀ w ∈ T, ∀ w' ∈ T, σ w = σ w' → w = w')
    (hsurj : ∀ v ∈ T, ∃ w, ∃ _ : w ∈ T, σ w = v)
    (hshift : ∀ w ∈ T, w + w₀ - σ w ∈ P.lattice) (z : ℂ) :
    veluAnalyticX P T c (z + w₀) = veluAnalyticX P T c z := by
  simp only [veluAnalyticX]
  congr 1
  rw [← veluAnalyticSum_shift P T w₀ σ hmem hinj hsurj hshift z]
  exact Finset.sum_congr rfl fun w _ => by ring_nf

/-- ★★★★★★★★★★★★★★★★★★★★★★★★**Vélu の公式の解析側——Liouville で閉じる形**。

    `℘_{Λ′}(z) = Σ_{w ∈ T} ℘_Λ(z + w) − c`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★仮定は 3 つ:

* `hper`  : 右辺が `Λ′`-周期的（★`veluAnalyticX_shift` で取れる——**代表系の置換だけ**）
* `hdiff` : ☆**差が整である**（極が打ち消し合うこと——極の解析が残る）
* `h0`    : ☆差が原点で `0`

★★★`hper` は本ファイルで無条件に取れており、残るのは `hdiff`（極の打ち消し）と `h0`。
☆すなわち **Vélu の公式の解析側は「極の解析」1 点に絞られた**。

★★★★これが `Found/GenEll/Velu.lean`（第 586-593）の代数側と対になる。
両者を繋げば `Found/GaloisRep/VeluNormalized.lean` の
`htFalt_isogeny_le_of_analytic`（第 596）の入力が揃う。 -/
theorem weierstrassP_eq_of_liouville (P P' : PeriodPair) (T : Finset ℂ) (c : ℂ)
    (hper : ∀ z : ℂ, ∀ l ∈ P'.lattice, veluAnalyticX P T c (z + l) = veluAnalyticX P T c z)
    (hdiff : Differentiable ℂ (fun z => P'.weierstrassP z - veluAnalyticX P T c z))
    (h0 : P'.weierstrassP 0 - veluAnalyticX P T c 0 = 0) (z : ℂ) :
    P'.weierstrassP z = veluAnalyticX P T c z := by
  have hkey : ∀ (y : ℂ), ∀ l ∈ P'.lattice,
      (fun z => P'.weierstrassP z - veluAnalyticX P T c z) (y + l)
        = (fun z => P'.weierstrassP z - veluAnalyticX P T c z) y := by
    intro y l hl
    simp only
    rw [P'.weierstrassP_add_coe y ⟨l, hl⟩, hper y l hl]
  have h := elliptic_liouville_eq_zero P' _ hdiff hkey h0 z
  simpa [sub_eq_zero] using h

def veluAnalyticX_shift.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の X は Λ′-周期的——代表系は平行移動で置換される。★無条件)",
    sectionId := "genell-lemma-3-5" }

def weierstrassP_eq_of_liouville.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の公式の解析側——Liouville で閉じる形。残るのは極の解析)",
    sectionId := "genell-lemma-3-5" }

def weierstrassP_eq_of_liouville.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆仮定 hdiff: ℘_{Λ′}(z) − Σ_{w∈T} ℘_Λ(z+w) + c が整であること。" ++
       "★両辺とも Λ′ の各点で 2 位の極をもち、主要部が一致するので差は整になる。" ++
       "★★mathlib は order_weierstrassP(各格子点で 2 位の極)を持っているので" ++
       "道具はある(2026-08-29 に測定)") 8,
    .implicitStep
      ("☆仮定 h0: 差が原点で 0 であること。★c の取り方(c = Σ_{w≠0} ℘_Λ(w))で決まる") 6,
    .implicitStep
      ("★★★到達点(2026-08-29、第 599): Vélu の公式の解析側が" ++
       "「代表系の置換」(無条件で取れた)と「極の解析」(残り)に分離した。" ++
       "★代数側は Found/GenEll/Velu.lean が持っている(第 586-593)") 8 ]

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★極が打ち消し合うこと -/

/-- ★★★★★★★★★★**`℘` から主要部を引くと格子点で解析的**。

    `℘_Λ(z) − 1/(z−p)²` は `p ∈ Λ` で解析的

★mathlib の `weierstrassPExcept`（`l₀` 項を抜いた `℘`）がそのまま使える:
`℘[Λ] z = ℘[Λ − p] z + (1/(z−p)² − 1/p²)`。 -/
theorem weierstrassP_sub_pole_analyticAt (P : PeriodPair) (p : ℂ) (hp : p ∈ P.lattice) :
    AnalyticAt ℂ (fun z => P.weierstrassP z - 1 / (z - p) ^ 2) p := by
  have hEq : (fun z => P.weierstrassP z - 1 / (z - p) ^ 2)
      = fun z => P.weierstrassPExcept p z - 1 / p ^ 2 := by
    funext z; rw [← P.weierstrassPExcept_add ⟨p, hp⟩]; ring
  rw [hEq]
  exact (((P.differentiableOn_weierstrassPExcept p).analyticOnNhd
    P.isOpen_compl_lattice_sdiff) p (by simp)).sub analyticAt_const

/-- ★★★★★★★★平行移動した `℘` からも同じ主要部を引ける——`(z+w) − (p+w) = z − p` だから。 -/
theorem shifted_sub_pole_analyticAt (P : PeriodPair) (p w : ℂ) (h : p + w ∈ P.lattice) :
    AnalyticAt ℂ (fun z => P.weierstrassP (z + w) - 1 / (z - p) ^ 2) p := by
  have hf : AnalyticAt ℂ (fun z : ℂ => z + w) p := analyticAt_id.add analyticAt_const
  have hg : AnalyticAt ℂ (fun z => P.weierstrassP z - 1 / (z - (p + w)) ^ 2)
      ((fun z : ℂ => z + w) p) := weierstrassP_sub_pole_analyticAt P (p + w) h
  refine (AnalyticAt.comp (f := fun z : ℂ => z + w) (x := p) hg hf).congr ?_
  filter_upwards with z
  simp only [Function.comp_apply]
  ring_nf

/-- ★★★★★格子の外では平行移動した `℘` は解析的。 -/
theorem shifted_analyticAt (P : PeriodPair) (p w : ℂ) (h : p + w ∉ P.lattice) :
    AnalyticAt ℂ (fun z => P.weierstrassP (z + w)) p := by
  have hf : AnalyticAt ℂ (fun z : ℂ => z + w) p := analyticAt_id.add analyticAt_const
  have hg : AnalyticAt ℂ P.weierstrassP ((fun z : ℂ => z + w) p) :=
    P.analyticOnNhd_weierstrassP (p + w) h
  exact AnalyticAt.comp (f := fun z : ℂ => z + w) (x := p) hg hf

/-- ★★★★★★★★★★★★★★★★★★**極が打ち消し合う**——`p ∈ Λ′` のとき。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`Λ′` の点 `p` では、`℘_{Λ′}` と（代表系のうちちょうど 1 つの）`℘_Λ(·+w₀)` が
**同じ主要部 `1/(z−p)²`** を持つので、差は解析的になる。 -/
theorem veluDiff_analyticAt_of_mem (P P' : PeriodPair) (T : Finset ℂ) (c p : ℂ)
    (hp : p ∈ P'.lattice) (w₀ : ℂ) (hw₀ : w₀ ∈ T) (hpw₀ : p + w₀ ∈ P.lattice)
    (hother : ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice) :
    AnalyticAt ℂ (fun z => P'.weierstrassP z - veluAnalyticX P T c z) p := by
  have hEq : (fun z => P'.weierstrassP z - veluAnalyticX P T c z)
      = fun z => ((P'.weierstrassP z - 1 / (z - p) ^ 2)
            - (P.weierstrassP (z + w₀) - 1 / (z - p) ^ 2))
          - (∑ w ∈ T.erase w₀, P.weierstrassP (z + w)) + c := by
    funext z
    simp only [veluAnalyticX]
    rw [← Finset.add_sum_erase T (fun w => P.weierstrassP (z + w)) hw₀]
    ring
  rw [hEq]
  refine (((weierstrassP_sub_pole_analyticAt P' p hp).sub
    (shifted_sub_pole_analyticAt P p w₀ hpw₀)).sub ?_).add analyticAt_const
  refine Finset.analyticAt_fun_sum _ fun w hw => ?_
  exact shifted_analyticAt P p w
    (hother w (Finset.mem_of_mem_erase hw) (Finset.ne_of_mem_erase hw))

/-- ★★★★★★★★★★★★格子の外では両方とも解析的。 -/
theorem veluDiff_analyticAt_of_notMem (P P' : PeriodPair) (T : Finset ℂ) (c p : ℂ)
    (hp : p ∉ P'.lattice) (hall : ∀ w ∈ T, p + w ∉ P.lattice) :
    AnalyticAt ℂ (fun z => P'.weierstrassP z - veluAnalyticX P T c z) p := by
  have hEq : (fun z => P'.weierstrassP z - veluAnalyticX P T c z)
      = fun z => (P'.weierstrassP z - ∑ w ∈ T, P.weierstrassP (z + w)) + c := by
    funext z; simp only [veluAnalyticX]; ring
  rw [hEq]
  refine (((P'.analyticOnNhd_weierstrassP p hp).sub ?_)).add analyticAt_const
  exact Finset.analyticAt_fun_sum _ fun w hw => shifted_analyticAt P p w (hall w hw)

/-- ★★★★★★`p ∉ Λ′` なら `p + w ∉ Λ`（`w ∈ T ⊆ Λ′`・`Λ ⊆ Λ′` だから）。 -/
theorem notMem_of_notMem (P P' : PeriodPair) (hle : P.lattice ≤ P'.lattice) (T : Finset ℂ)
    (hT : ∀ w ∈ T, w ∈ P'.lattice) (p : ℂ) (hp : p ∉ P'.lattice) (w : ℂ) (hw : w ∈ T) :
    p + w ∉ P.lattice := by
  intro hc
  refine hp ?_
  have h1 : (p + w) - w ∈ P'.lattice := P'.lattice.sub_mem (hle hc) (hT w hw)
  simpa using h1

/-- ★★★★★★★★★★★★★★★★★★★★★★**差は整である**——★極の打ち消しが済んだ形。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★仮定は **`T` が `Λ′/Λ` の代表系である**ことだけ（`hT`・`hrep`）:

* `hT`  : `T ⊆ Λ′`
* `hrep`: `p ∈ Λ′` ならちょうど 1 つの `w₀ ∈ T` が `p + w₀ ∈ Λ` を満たす

☆すなわち `weierstrassP_eq_of_liouville`（第 599）の仮定 `hdiff` は
**代表系の性質だけから従う**。★残るのは `h0`（差が原点で `0`）＝`c` の取り方である。 -/
theorem veluDiff_differentiable (P P' : PeriodPair) (hle : P.lattice ≤ P'.lattice)
    (T : Finset ℂ) (c : ℂ) (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice) :
    Differentiable ℂ (fun z => P'.weierstrassP z - veluAnalyticX P T c z) := by
  intro p
  by_cases hp : p ∈ P'.lattice
  · obtain ⟨w₀, hw₀, hpw₀, hother⟩ := hrep p hp
    exact (veluDiff_analyticAt_of_mem P P' T c p hp w₀ hw₀ hpw₀ hother).differentiableAt
  · exact (veluDiff_analyticAt_of_notMem P P' T c p hp
      (fun w hw => notMem_of_notMem P P' hle T hT p hp w hw)).differentiableAt

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★**Vélu の公式の解析側——代表系だけから**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

    `℘_{Λ′}(z) = Σ_{w ∈ T} ℘_Λ(z + w) − c`

★★★仮定は **`T` が `Λ′/Λ` の代表系であること**（`hT`・`hrep`）、
**平行移動が代表系を置換すること**（`hshift`）、
そして `h0`（`c` の取り方で決まる正規化）だけになった。

☆☆**極の解析は済んだ**（`veluDiff_differentiable`）——第 599 で残っていた `hdiff` である。 -/
theorem weierstrassP_eq_veluAnalyticX (P P' : PeriodPair) (hle : P.lattice ≤ P'.lattice)
    (T : Finset ℂ) (c : ℂ) (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    (hper : ∀ z : ℂ, ∀ l ∈ P'.lattice, veluAnalyticX P T c (z + l) = veluAnalyticX P T c z)
    (h0 : P'.weierstrassP 0 - veluAnalyticX P T c 0 = 0) (z : ℂ) :
    P'.weierstrassP z = veluAnalyticX P T c z :=
  weierstrassP_eq_of_liouville P P' T c hper
    (veluDiff_differentiable P P' hle T c hT hrep) h0 z

def veluDiff_differentiable.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(極が打ち消し合う——代表系の性質だけから従う。★無条件)",
    sectionId := "genell-lemma-3-5" }

def weierstrassP_eq_veluAnalyticX.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の公式の解析側——代表系と正規化だけから)",
    sectionId := "genell-lemma-3-5" }

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★正規化定数と最終形 -/

/-- ★★★★★**`℘(0) = 0`**——mathlib の定義（`∑' l, (1/(z−l)² − 1/l²)`）では
`z = 0` の各項が `1/l² − 1/l² = 0` になる（`l = 0` の項も junk value で `0`）。

★これは「極での値」ではなく**除去可能特異点を埋めた関数の値**として整合している
（`weierstrassPExcept_add` が `z = l₀` でも成り立つのはそのためである）。 -/
theorem weierstrassP_zero (P : PeriodPair) : P.weierstrassP 0 = 0 := by
  simp [PeriodPair.weierstrassP]

/-- ★★★★★★**Vélu の正規化定数** `c = Σ_{w ∈ T∖{0}} ℘_Λ(w)`。

★原文の形 `℘(z) + Σ_{w≠0}[℘(z+w) − ℘(w)]` の第 2 項の定数部分である。 -/
noncomputable def veluAnalyticC (P : PeriodPair) (T : Finset ℂ) : ℂ :=
  ∑ w ∈ T.erase 0, P.weierstrassP w

/-- ★★★★★★★★**正規化定数を入れると原点で `0`**。 -/
theorem veluAnalyticX_zero (P : PeriodPair) (T : Finset ℂ) (h0T : (0:ℂ) ∈ T) :
    veluAnalyticX P T (veluAnalyticC P T) 0 = 0 := by
  simp only [veluAnalyticX, veluAnalyticC]
  rw [← Finset.add_sum_erase T (fun w => P.weierstrassP (0 + w)) h0T]
  simp

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**Vélu の公式の解析側（最終形）**

    `℘_{Λ′}(z) = Σ_{w ∈ T} ℘_Λ(z + w) − Σ_{w ∈ T∖{0}} ℘_Λ(w)`

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★**仮定は `T` が `Λ′/Λ` の代表系であることだけ**である:

| 仮定 | 内容 |
|---|---|
| `hle` | `Λ ⊆ Λ′` |
| `h0T` | `0 ∈ T`（自明な剰余類の代表） |
| `hT` | `T ⊆ Λ′` |
| `hrep` | `p ∈ Λ′` ならちょうど 1 つの `w₀ ∈ T` が `p + w₀ ∈ Λ` |
| `hper` | 平行移動が代表系を置換する（★`veluAnalyticX_shift` で取れる） |

★★★★**極の解析も正規化も済んだ**——第 598（Liouville）・第 599（周期性）・
第 600（極の打ち消し）・第 601（正規化定数）で塞いだ。

☆残るのは `hper` を代表系の定義から出すこと（`veluAnalyticX_shift` に `σ` を与えること）と、
この `℘` の等式を `Found/GenEll/Velu.lean` の**代数側の `veluXGen`** に翻訳することである
——後者には `℘` の加法定理（mathlib に無い）が要る。 -/
theorem weierstrassP_eq_velu (P P' : PeriodPair) (hle : P.lattice ≤ P'.lattice)
    (T : Finset ℂ) (h0T : (0:ℂ) ∈ T) (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    (hper : ∀ z : ℂ, ∀ l ∈ P'.lattice,
      veluAnalyticX P T (veluAnalyticC P T) (z + l)
        = veluAnalyticX P T (veluAnalyticC P T) z)
    (z : ℂ) :
    P'.weierstrassP z = veluAnalyticX P T (veluAnalyticC P T) z :=
  weierstrassP_eq_veluAnalyticX P P' hle T _ hT hrep hper
    (by rw [weierstrassP_zero, veluAnalyticX_zero P T h0T]; ring) z

def weierstrassP_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(℘(0) = 0——除去可能特異点を埋めた値。★無条件)",
    sectionId := "genell-lemma-3-5" }

def weierstrassP_eq_velu.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の公式の解析側・最終形——仮定は代表系であることだけ)",
    sectionId := "genell-lemma-3-5" }

def weierstrassP_eq_velu.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆仮定 hper は veluAnalyticX_shift(第 599)で取れるが、" ++
       "代表系の定義から σ(平行移動が誘導する置換)を作る段が残る") 5,
    .implicitStep
      ("☆この ℘ の等式を Found/GenEll/Velu.lean の代数側 veluXGen(第 591)に" ++
       "翻訳すること。★℘ の加法定理が要る(mathlib に無い、2026-08-29 に測定)") 8,
    .implicitStep
      ("★★★到達点(2026-08-29、第 601): Vélu の公式の解析側が" ++
       "「T が Λ′/Λ の代表系である」だけから従う形になった。" ++
       "★極の解析(第 600)・Liouville(第 598)・周期性(第 599)・正規化(第 601)で塞いだ") 9 ]

end ABC3.Found.GenEll
