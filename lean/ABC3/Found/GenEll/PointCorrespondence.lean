/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.GlobalExtend
import ABC3.Found.GenEll.Remark151Sigma

/-!
# [GenEll] Remark 1.5.1 —— **点の対応 `ePt` を作る**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

## ★★★★★★★★★★仮定が定理になった

`Remark151Sigma.lean` の `remark_1_5_1_bdeq` は
**点の対応 `ePt` とその両立 `hcompat` をデータとして受けて**いた。
★`GlobalExtend.lean` で `X′(F) ≃ X′(𝓞_F)`（`ℤ` 上固有な `X′`）が取れたので、
**`ePt` は構成できる**。

## ★★★★★★機構は 2 段

1. `F`-点を作る —— `Spec F ⟶ Spec 𝓞_F[1/n!] ⟶ X_n --φ--> X′_n ⟶ X′`。
2. 固有性で `Spec 𝓞_F` へ延ばす（`exists_global_extend`）。

★★両立 `hcompat` は**生成点で一致すること**から出る:
`Spec 𝓞_F[1/n!]` は整、`X′` は分離的、生成点は稠密
——`ext_of_isDominant_of_isSeparated` 1 本。

★★★`Spec 𝓞_F[1/n!] ⟶ Spec 𝓞_F` の分数体が `F` であること（`IsFractionRing A F`）が
稠密性の入力である。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Limits IsDedekindDomain

variable (F : Type) [Field F] [NumberField F]

/-! ## ★★★★★★★★★★点の対応を構成する -/

/-- ★★★★★★★★★★**点の対応 `ePt` を構成する**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

段 `n` での同型 `φ` から、`Spec 𝓞_F`-点の対応が**作れる**。

★★機構: `F`-点を作って固有性で延ばす。両立は生成点で一致することから
（`ext_of_isDominant_of_isSeparated`）。 -/
theorem exists_ePt {X X' : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ)) [IsProper f']
    {n : ℕᵒᵖ}
    (A : Type) [CommRing A] [IsDomain A] [Algebra (NumberField.RingOfIntegers F) A]
    [Algebra A F] [IsScalarTower (NumberField.RingOfIntegers F) A F] [IsFractionRing A F]
    (hinv : IsUnit (algebraMap ℤ A ((Nat.factorial n.unop : ℕ) : ℤ)))
    (φ : bcObj f n ⟶ bcObj f' n) :
    ∃ ePt : (specRingOfIntegers F ⟶ X) → (specRingOfIntegers F ⟶ X'),
      ∀ xF, Spec.map (CommRingCat.ofHom (algebraMap (NumberField.RingOfIntegers F) A)) ≫ ePt xF
        = liftPointToBc F A hinv f xF ≫ φ ≫
          pullback.snd (overRatTowerDiagram.obj n).hom f' := by
  classical
  choose ePt hePt using fun xF : specRingOfIntegers F ⟶ X =>
    exists_global_extend F f' (Spec.map (CommRingCat.ofHom (algebraMap A F)) ≫
      (liftPointToBc F A hinv f xF ≫ φ ≫ pullback.snd (overRatTowerDiagram.obj n).hom f'))
  refine ⟨ePt, fun xF => ?_⟩
  haveI : IsDominant (Spec.map (CommRingCat.ofHom (algebraMap A F))) :=
    isDominant_specMap_fractionRing A F
  refine ext_of_isDominant_of_isSeparated f' (specZIsTerminal.hom_ext _ _)
    (Spec.map (CommRingCat.ofHom (algebraMap A F))) ?_
  rw [← Category.assoc, ← Spec.map_comp, ← CommRingCat.ofHom_comp,
    ← IsScalarTower.algebraMap_eq]
  exact hePt xF

/-! ## ★★★★★★★★★★到達点 —— `Remark 1.5.1` の後半（`ePt` は構成） -/

/-- ★★★★★★★★★★**[GenEll] Remark 1.5.1 の後半 —— `ePt` を構成した形**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

`DivisorDescent.lean` の `exists_pair_descent` が出す
「有限段 `n` での同型＋可換正方形」から、
**点の対応 `ePt` を作り**、`log-cond` の BD-同値を出す。

★★★**もう `ePt` は仮定ではない**——固有性から構成される
（`GlobalExtend.lean` の `exists_global_extend`）。
★★`extend_unique`（`DedekindExtend.lean`）により、その `ePt` は**一意**でもある。

★定数は `∑_{q ∣ n!} log q`——原文の `∑_{q ∈ Σ} log q` そのもの。

★★★★残る仮定 `hI` / `hI′` は「点が因子を通らない」ことであり、
原文が定義域を `U_X(ℚ̄) = X ∖ D` としていることに対応する。 -/
theorem remark_1_5_1_of_descent {Z Z' X X' : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ))
    (iZ : Z ⟶ X) (iZ' : Z' ⟶ X') [IsClosedImmersion iZ] [IsClosedImmersion iZ']
    [IsProper f']
    {n : ℕᵒᵖ}
    (A : Type) [CommRing A] [IsDomain A] [Algebra (NumberField.RingOfIntegers F) A]
    [IsLocalization
      (Submonoid.powers ((Nat.factorial n.unop : ℕ) : NumberField.RingOfIntegers F)) A]
    [Algebra A F] [IsScalarTower (NumberField.RingOfIntegers F) A F] [IsFractionRing A F]
    (hinv : IsUnit (algebraMap ℤ A ((Nat.factorial n.unop : ℕ) : ℤ)))
    (φ : bcObj f n ⟶ bcObj f' n) (ψ : bcObj (iZ ≫ f) n ⟶ bcObj (iZ' ≫ f') n)
    [IsIso φ] [IsIso ψ]
    (hsq : ψ ≫ bcBC (iZ' ≫ f') f' iZ' n = bcBC (iZ ≫ f) f iZ n ≫ φ)
    (ch : FinitePlace F → ℕ) (hchprime : ∀ v, (ch v).Prime)
    (hover : ∀ v : FinitePlace F,
      (v.asIdeal).LiesOver (Ideal.span {((ch v : ℕ) : ℤ)}))
    (hI : ∀ xF, pullbackIdeal F iZ.ker xF ≠ 0) :
    ∃ ePt : (specRingOfIntegers F ⟶ X) → (specRingOfIntegers F ⟶ X'),
      (∀ xF, pullbackIdeal F iZ'.ker (ePt xF) ≠ 0) →
      BDeq (fun xF => logCond F iZ.ker xF) (fun xF => logCond F iZ'.ker (ePt xF)) := by
  obtain ⟨ePt, hcompat⟩ := exists_ePt F f f' A hinv φ
  exact ⟨ePt, fun hI' => remark_1_5_1_bdeq F f f' iZ iZ' A hinv φ ψ hsq ch hchprime hover
    ePt hcompat hI hI'⟩

/-! ### ★出典の紐付け(`.src`)

★★**まだ条つきである。** 残っているのは「`ℚ`-同型から降下データ `(n, φ, ψ, hsq)` を取る」段
——`DivisorDescent.lean` の `exists_pair_descent` がそれを出すので、
**あとは呼ぶだけ**である（本ファイルではまだ呼んでいない）。 -/

def exists_ePt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(点の対応 ePt を固有性から構成する)",
    sectionId := "genell-rem-1-5-1" }

def remark_1_5_1_of_descent.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(後半——降下データを与えられた形。ℚ-同型から取る段は含まない)",
    sectionId := "genell-rem-1-5-1" }

def remark_1_5_1_of_descent.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_global_extend(固有性から点が延びる)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_global_extend") 8,
    .citation "[ABC3]" "remark_1_5_1_bdeq(ePt を受けた形の BD-同値)"
      (.inProject "ABC3" "ABC3.Found.GenEll.remark_1_5_1_bdeq") 9,
    .citation "[ABC3]" "exists_pair_descent(ℚ-同型から降下データ)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_pair_descent") 9,
    .implicitStep
      ("★★★★★点の対応 ePt は**もう仮定ではない**——固有性から構成される。" ++
       "★extend_unique により一意でもある") 9,
    .implicitStep
      ("★★残る仮定 hI / hI′ は「点が因子を通らない」ことであり、" ++
       "原文が定義域を U_X(ℚ̄) = X ∖ D としていることに対応する") 9,
    .implicitStep
      ("★★★残る段: ℚ-同型から降下データ (n, φ, ψ, hsq) を取る。" ++
       "exists_pair_descent がそれを出すので、あとは呼ぶだけである") 9 ]

end ABC3.Found.GenEll
