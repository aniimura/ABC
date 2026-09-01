/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateLevel
import ABC3.Found.GaloisRep.TateWiring
import ABC3.Meta.Claim

/-!
# 第 1203 ブロック —— **`tateProj` の核はちょうど `l^n · T_l E`**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★これは何か——`T_l E / l ≅ E[l]`

第 1201-1202 で「第 1 層の射影は位数 `l` の点を与え、逆に `l`-捩れ点は
第 1 層に持ち上がる」を取った。本ブロックは**核**を取る:

    (tateProj W l n f : Point) = 0  ↔  ∃ g, (l^n) • g = f

★★★これで `T_l E / l^n · T_l E ≅ E[l^n]` が言える。
`n = 1` が `HasLCyclicJ`（`mod l` の直線）を
**点の直線**へ一意に対応させる段である。

## ★★★★★★★★測ったこと——`tateModule` と `limTors` は**定義から同じ**

`Found/GaloisRep/TateLevel.lean` の `exists_smul_of_proj_zero` は
抽象群 `A` の `limTors A l`（`(nsmulHom A (l^m)).ker` の逆極限）で書かれており、
`Interface/GaloisRep/Torsion.lean` の `tateModule W l` は
`torsionPoints W (l^m)` の逆極限である。

☆両者の membership は `mem_ker_nsmulHom` が `Iff.rfl` なので**定義から等しく**、
`exact exists_smul_of_proj_zero l n f h` が**そのまま通る**（2026-09-02 実測）。
★したがって片方向は在庫の使い回しであり、新しいのは逆向きと `iff` の形である。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep ABC3.Meta

variable {F : Type} [Field F] [DecidableEq F]

/-! ## ★★★★★★★★★★核に入るなら `l^n` 倍で書ける -/

/-- ★★★★★★★★★★**第 `n` 射影の核は `l^n · T_l E` に入る**——★**無条件**（第 1203）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`tateModule W l` は `limTors W.toAffine.Point l` と**定義から同じ**なので、
在庫の `exists_smul_of_proj_zero` がそのまま当たる。 -/
theorem exists_smul_of_tateProj_eq_zero (W : WeierstrassCurve F) (l n : ℕ)
    (f : tateModule W l) (h : ((tateProj W l n f : W.toAffine.Point)) = 0) :
    ∃ g : tateModule W l, (l ^ n) • g = f :=
  exists_smul_of_proj_zero l n f h

/-! ## ★★★★★★逆向き -/

/-- ★★★★★★**`l^n` 倍は第 `n` 射影で消える**——★**無条件**（第 1203）。

☆`tateProj W l n g` は `E[l^n]` に入っているからである。 -/
theorem tateProj_smul_pow_eq_zero (W : WeierstrassCurve F) (l n : ℕ)
    (g : tateModule W l) :
    ((tateProj W l n ((l ^ n) • g) : W.toAffine.Point)) = 0 := by
  have h : tateProj W l n ((l ^ n) • g) = (l ^ n) • tateProj W l n g := map_nsmul _ _ _
  rw [h]
  simpa using (tateProj W l n g).2

/-! ## ★★★★★★★★★★★★★★★★`iff` の形 -/

/-- ★★★★★★★★★★★★★★★★
**`T_l E / l^n · T_l E ≅ E[l^n]`（核の記述）**——★**無条件**（第 1203）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★`n = 1` が `HasLCyclicJ`（`mod l` の直線）を**点の直線**へ
一意に対応させる段である。 -/
theorem tateProj_eq_zero_iff (W : WeierstrassCurve F) (l n : ℕ) (f : tateModule W l) :
    ((tateProj W l n f : W.toAffine.Point)) = 0 ↔ ∃ g : tateModule W l, (l ^ n) • g = f := by
  constructor
  · exact exists_smul_of_tateProj_eq_zero W l n f
  · rintro ⟨g, rfl⟩
    exact tateProj_smul_pow_eq_zero W l n g

/-- ★★★★★★★★★★★★**第 1 層の形**——★**無条件**（第 1203）。

☆`tateProj W l 1 f = 0 ↔ f ∈ l · T_l E`。
★★★これが「`mod l` で `≠ 0`」と「第 1 層の点が `≠ 0`」を結ぶ。 -/
theorem tateProj_one_eq_zero_iff (W : WeierstrassCurve F) (l : ℕ) (f : tateModule W l) :
    ((tateProj W l 1 f : W.toAffine.Point)) = 0 ↔ ∃ g : tateModule W l, l • g = f := by
  rw [tateProj_eq_zero_iff]
  simp

/-! ## ★出典の紐付け(`.src`) -/

def exists_smul_of_tateProj_eq_zero.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(第 n 射影の核は l^n · T_l E に入る。★無条件)",
    sectionId := "genell-thm-3-8" }

def tateProj_smul_pow_eq_zero.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(l^n 倍は第 n 射影で消える。★無条件)",
    sectionId := "genell-thm-3-8" }

def tateProj_eq_zero_iff.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(T_l E / l^n · T_l E ≅ E[l^n]——核の記述。★無条件)",
    sectionId := "genell-thm-3-8" }

def tateProj_eq_zero_iff.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_smul_of_proj_zero(TateLevel.lean、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_smul_of_proj_zero") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1203）**——`tateModule W l` は" ++
       "`limTors W.toAffine.Point l` と**定義から同じ**である" ++
       "（`mem_ker_nsmulHom` が `Iff.rfl`、2026-09-02 実測）ので、" ++
       "片方向は在庫がそのまま当たる。☆新しいのは逆向きと `iff` の形である。" ++
       "★`n = 1` が `HasLCyclicJ`（`mod l` の直線）を**点の直線**へ" ++
       "一意に対応させる段である。") 2 ]

def tateProj_one_eq_zero_iff.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(第 1 層の射影の核は l · T_l E。★無条件)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
