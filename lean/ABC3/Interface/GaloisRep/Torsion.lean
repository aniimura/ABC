import ABC3.Meta.Claim
import Mathlib.AlgebraicGeometry.EllipticCurve.Affine.Point
import Mathlib.FieldTheory.IsAlgClosed.Basic
import Mathlib.NumberTheory.Padics.PadicIntegers

/-!
# Galois 表現のスケルトン(1/3)—— **捩れ点と Tate 加群**

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★なぜ割るのか

`Interface/GenEll/GaloisRep.lean` の `TorsionGaloisRepData` は
`Curve : Type` / `EllClass : Type` と**世界ごと抽象化**していた。
★それでは「何を作れば埋まるか」が個数として見えない。

★★本ファイル群は、mathlib の**実在する楕円曲線**
(`WeierstrassCurve` + `IsElliptic`)の上で obligation を並べ直す。
★★★**1 件ずつ独立に埋められる。**

## ★★★捩れ部分群は posit しない

`E[n]` そのものは **mathlib の点の群から今すぐ作れる**——
`W.toAffine.Point` は `AddCommGroup` を持つからである。
★したがって本ファイルは `torsionPoints` を **`def` として定義し**、
**構造定理だけ**を `Interface` で受ける。

## ★★本ファイルが受ける 2 件

| # | 受けるもの | mathlib / FLT の在庫(2026-08-17 実測) |
|---|---|---|
| G1 | `E[n] ≅ (ℤ/n)²`(代数閉体上、`n` が可逆) | ★**無い**。FLT の `Torsion.lean` にもあるが `sorry` |
| G2 | Tate 加群 `T_l E ≅ ℤ_l²` | ★**無い**(`TateModule` は mathlib に 0 件) |
-/

namespace ABC3.Interface.GaloisRep

open ABC3.Meta WeierstrassCurve

/-! ## ★★捩れ部分群(これは posit ではない) -/

/-- **`E[n]`** —— `n` 倍で消える点。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★**これは今すぐ書ける**——mathlib の `W.toAffine.Point` が
`AddCommGroup` を持つからである。★posit するのは**構造定理の方だけ**である。 -/
def torsionPoints {K : Type} [Field K] [DecidableEq K] (W : WeierstrassCurve K) (n : ℕ) :
    AddSubgroup W.toAffine.Point where
  carrier := {P : W.toAffine.Point | n • P = 0}
  add_mem' := by
    intro a b ha hb
    simp only [Set.mem_setOf_eq] at *
    rw [smul_add, ha, hb, add_zero]
  zero_mem' := by simp
  neg_mem' := by
    intro a ha
    simp only [Set.mem_setOf_eq] at *
    rw [smul_neg, ha, neg_zero]

/-! ## ★★★G1 —— `E[n] ≅ (ℤ/n)²` -/

/-- **(G1)** 代数閉体上の `n`-捩れの構造 `E[n] ≅ (ℤ/n)²`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★**`GL₂` が出るのはこれがあるからである。**捩れが階数 2 だから
自己同型群が `GL₂` になる。★これが無いと表現の**行き先が書けない**。

★`n` は体の標数と素であることが要る(そうでないと階数が落ちる)。 -/
structure TorsionStructureData where
  /-- ★★**代数閉体上、`n` が可逆なら `(ℤ/n)²` と同型**。 -/
  structure_eq : ∀ (K : Type) [Field K] [DecidableEq K] [IsAlgClosed K] (W : WeierstrassCurve K),
    W.IsElliptic → ∀ n : ℕ, 0 < n → (∀ k : ℕ, 1 ≤ k → k ≤ n → (k : K) ≠ 0) →
    Nonempty (torsionPoints W n ≃+ (ZMod n × ZMod n))
  /-- ★`n` 捩れは有限。 -/
  torsion_finite : ∀ (K : Type) [Field K] [DecidableEq K] [IsAlgClosed K] (W : WeierstrassCurve K),
    W.IsElliptic → ∀ n : ℕ, 0 < n → (n : K) ≠ 0 → Finite (torsionPoints W n)

def TorsionStructureData.waiting : WaitingFor :=
  { what := "(G1) 代数閉体上の n-捩れの構造定理 E[n] ≅ (Z/n)^2(n が標数と素のとき)"
    trackB := "Found/GaloisRep — ★mathlib は `WeierstrassCurve.Affine.Point` の `AddCommGroup` と分点多項式(`DivisionPolynomial`)までを持つが、**構造定理は無い**(2026-08-16 実測)。★★**FLT プロジェクトにも無い**——`FLT/EllipticCurve/Torsion.lean` は 124 行中 sorry 10 件で `n_torsion_card` / `n_torsion_dimension` が**いずれも sorry**(2026-08-17、clone して計数)。★★★これが S3 の最初の壁である" }

/-! ## ★★G2 —— Tate 加群 -/

/-- **`T_l E = lim_n E[l^n]`** —— `l` 進 Tate 加群。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★これも posit ではない

逆極限は `E[l^n]` の直積の**部分群**として今すぐ書ける——
「`l` 倍で一つ下の層に落ちる列」の全体である。
★★posit するのは **`ℤ_l²` と同型であること**だけである。

★★★§9-401 の教訓: `tateModule` を素の `Type` として posit し
`proj` に条件を課さないと、`tateModule := ℤ_l × ℤ_l`、`proj := 0` で
**空虚に埋まってしまう**。★定義できるものは定義する。 -/
def tateModule {K : Type} [Field K] [DecidableEq K] (W : WeierstrassCurve K) (l : ℕ) :
    AddSubgroup (∀ n : ℕ, torsionPoints W (l ^ n)) where
  carrier := {f | ∀ n : ℕ, l • ((f (n + 1) : W.toAffine.Point)) = (f n : W.toAffine.Point)}
  add_mem' := by
    intro a b ha hb n
    simp only [Set.mem_setOf_eq, Pi.add_apply, AddSubgroup.coe_add, smul_add] at *
    rw [ha n, hb n]
  zero_mem' := by
    intro n
    simp
  neg_mem' := by
    intro a ha n
    simp only [Set.mem_setOf_eq, Pi.neg_apply, AddSubgroup.coe_neg, smul_neg] at *
    rw [ha n]

/-- ★`T_l E → E[l^n]` の射影。 -/
def tateProj {K : Type} [Field K] [DecidableEq K] (W : WeierstrassCurve K) (l n : ℕ) :
    tateModule W l →+ torsionPoints W (l ^ n) where
  toFun := fun f => (f : ∀ m : ℕ, torsionPoints W (l ^ m)) n
  map_zero' := rfl
  map_add' := fun _ _ => rfl

/-- ★★射影は `l` 倍と可換である(逆極限であることの中身)。 -/
theorem tateProj_smul {K : Type} [Field K] [DecidableEq K] (W : WeierstrassCurve K) (l n : ℕ)
    (f : tateModule W l) :
    l • ((tateProj W l (n + 1) f : W.toAffine.Point)) = (tateProj W l n f : W.toAffine.Point) :=
  f.2 n

/-- **(G2)** `T_l E ≅ ℤ_l²`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★原文の `GL₂(ℤ_l)` の `ℤ_l²` がこれである。

★標数 0 を要求する: (G1) の個数勘定が `∀ k ≤ n, (k : K) ≠ 0` を要るので、
`n = l^m` を `m → ∞` で走らせるには標数 0 が要る(§9-398 の逸脱記録)。
★★ABC は `ℚ` 上(の代数閉包)で使うので十分である。 -/
structure TateModuleData where
  /-- 台となる捩れの構造定理。 -/
  toTorsionStructureData : TorsionStructureData
  /-- ★★**`T_l E` は階数 2 の自由 `ℤ_l` 加群**——**群同型**として。 -/
  freeRankTwo : ∀ (K : Type) [Field K] [DecidableEq K] [IsAlgClosed K] [CharZero K]
    (W : WeierstrassCurve K), W.IsElliptic → ∀ l : ℕ, ∀ _ : Fact l.Prime,
    Nonempty (tateModule W l ≃+ (ℤ_[l] × ℤ_[l]))

def TateModuleData.waiting : WaitingFor :=
  { what := "(G2) l 進 Tate 加群 T_l E = lim E[l^n] が階数 2 の自由 Z_l 加群であること"
    trackB := "Found/GaloisRep — ★mathlib に `TateModule` は **0 件**(2026-08-17 実測)。★(G1) に従属する——`E[l^n] ≅ (Z/l^n)^2` の**両立する基底**を取り、`PadicInt.lift` で逆極限を渡る" }

/-! ## ★出典の紐付け(`.src`) -/

def torsionPoints.src : Source :=
  { paper := "GenEll", pdfPage := 19, item := "Theorem 3.8(n-捩れ部分群の定義のみ)",
    sectionId := "genell-thm-3-8" }

def TorsionStructureData.src : Source :=
  { paper := "GenEll", pdfPage := 19, item := "Theorem 3.8(捩れの構造定理のみ——Galois 作用は含まない)",
    sectionId := "genell-thm-3-8" }

def TateModuleData.src : Source :=
  { paper := "GenEll", pdfPage := 19, item := "Theorem 3.8(Tate 加群のみ——Galois 作用は含まない)",
    sectionId := "genell-thm-3-8" }

end ABC3.Interface.GaloisRep
