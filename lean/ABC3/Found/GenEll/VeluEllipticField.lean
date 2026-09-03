/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluLatticePoint
import ABC3.Found.GenEll.VeluEllipticDescent
import ABC3.Found.GenEll.VeluPointSet
import ABC3.Found.GenEll.JSurjective
import ABC3.Meta.Claim

/-!
# 第 1334 ブロック —— **`ℂ` に埋め込める体の上で Vélu の商は楕円**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か——`VeluQuotOK` の楕円性の本体

第 1333 で `ℂ` の**格子曲線**については取れた。★本ブロックは

| 段 | 使うもの | 在庫 |
|---|---|---|
| 1 | `ℂ` へ送る（`Δ ≠ 0` は単射準同型で降りる） | `isElliptic_veluQuotientFull_of_map`（第 1331） |
| 2 | 点集合も `f` で送られる | `image_pointCoords_rhPoint_nsmul` |
| 3 | 一意化の変数変換 `C • (W⊗ℂ) = latticeCurve P` | `exists_periodPair_of_isElliptic` |
| 4 | 商は変数変換と可換 | `veluQuotientFull_vcPoint_eq`（第 974） |
| 5 | 格子曲線の商は楕円 | `isElliptic_veluQuotientFull_nsmul_lattice`（第 1333） |

を繋いで、**`ℂ` に埋め込める任意の体**の上で `⟨Q⟩` による Vélu の商が楕円であることを示す。

☆`hB`・`hBx`（点集合の反転安定性）は第 949 の在庫が自動で与える。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve ABC3.Meta

open scoped Classical

/-- ★★★★★★★★**変数変換した先が楕円なら元も楕円**——★**無条件**（第 1334）。

☆`C⁻¹ • (C • X) = X` である。 -/
theorem isElliptic_of_variableChange_smul {F : Type} [Field F] (C : VariableChange F)
    (X : WeierstrassCurve F) (h : (C • X).IsElliptic) : X.IsElliptic := by
  haveI := h
  have hx : C⁻¹ • (C • X) = X := inv_smul_smul C X
  rw [← hx]
  infer_instance

/-- ★★★★★★曲線が等しければ `⟨Q⟩` による商の楕円性は移る（第 1334）。 -/
theorem isElliptic_velu_congr_curve {F : Type} [Field F] {W₁ W₂ : WeierstrassCurve F}
    (h : W₁ = W₂) (Q : W₁.toAffine.Point) (l : ℕ) :
    (veluQuotientFull W₁
        (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).IsElliptic
      ↔ (veluQuotientFull W₂
        (((Finset.range l).erase 0).image
          (fun k : ℕ => pointCoords (k • (h ▸ Q : W₂.toAffine.Point))))).IsElliptic := by
  subst h
  exact Iff.rfl

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**`ℂ` に埋め込める体の上で `⟨Q⟩` による Vélu の商は楕円**——★**無条件**（第 1334）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これが `Skeleton/GenEll/AlphaBridge.lean` の `VeluQuotOK` が要求する
**楕円性そのもの**である（数体は `ℂ` に埋め込めるから）。 -/
theorem isElliptic_veluQuotientFull_nsmul_of_embed {F : Type} [Field F]
    (W : WeierstrassCurve F) [W.IsElliptic] (f : F →+* ℂ)
    {Q : W.toAffine.Point} {l : ℕ} (hl : 0 < l) (hQ : addOrderOf Q = l) :
    (veluQuotientFull W
      (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))).IsElliptic := by
  classical
  haveI hW1 : (W.map f).IsElliptic :=
    ⟨isUnit_iff_ne_zero.2 (by
      rw [WeierstrassCurve.map_Δ]
      exact (map_ne_zero_iff f f.injective).2 W.isUnit_Δ.ne_zero)⟩
  -- ★ 段 1・2: `ℂ` へ送る
  refine isElliptic_veluQuotientFull_of_map f W _ ?_
  rw [← image_pointCoords_rhPoint_nsmul f W hQ]
  have hQ₁ : addOrderOf (rhPoint f W Q) = l := by rw [addOrderOf_rhPoint, hQ]
  -- ★ 段 3: 一意化の変数変換
  obtain ⟨P, C, hCP⟩ := exists_periodPair_of_isElliptic (W.map f) hW1
  refine isElliptic_of_variableChange_smul C _ ?_
  -- ★ 段 4: 商は変数変換と可換
  rw [veluQuotientFull_vcPoint_eq C (W.map f) _ hQ₁ (by norm_num) rfl,
    isElliptic_velu_congr_curve hCP]
  -- ★ 段 5: 格子曲線の商は楕円（第 1333）
  refine isElliptic_veluQuotientFull_nsmul_lattice P (latticeDisc_ne_zero P) hl ?_
  rw [addOrderOf_congr_curve hCP, addOrderOf_vcPoint, hQ₁]

/-! ## ★出典の紐付け(`.src`) -/

def isElliptic_of_variableChange_smul.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(変数変換した先が楕円なら元も楕円。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isElliptic_veluQuotientFull_nsmul_of_embed.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(ℂ に埋め込める体の上で ⟨Q⟩ による Vélu の商は楕円。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isElliptic_veluQuotientFull_nsmul_of_embed.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isElliptic_veluQuotientFull_nsmul_lattice(第 1333、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isElliptic_veluQuotientFull_nsmul_lattice") 1,
    .citation "[ABC3]" "veluQuotientFull_vcPoint_eq(第 974、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.veluQuotientFull_vcPoint_eq") 1,
    .citation "[ABC3]" "exists_periodPair_of_isElliptic(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_periodPair_of_isElliptic") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1334）**——`Skeleton/GenEll/AlphaBridge.lean` の " ++
       "`VeluQuotOK` が要求する**楕円性そのもの**である。" ++
       "☆数体は `ℂ` に埋め込めるので、残るのは埋め込みを取る一行だけである。") 3 ]

end ABC3.Found.GenEll
