import ABC3.Interface.GaloisRep.Torsion
import Mathlib.LinearAlgebra.Matrix.GeneralLinearGroup.Defs
import Mathlib.FieldTheory.Galois.Basic
import Mathlib.NumberTheory.NumberField.Basic
import Mathlib.FieldTheory.IsAlgClosed.Basic
import Mathlib.NumberTheory.Padics.RingHoms

/-!
# Galois 表現のスケルトン(2/3)—— **表現そのものと像**

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★本ファイルが受ける 3 件

| # | 受けるもの | mathlib の在庫(2026-08-17 実測) |
|---|---|---|
| G3 | `ρ_{E,l} : Gal(K̄/K) → GL₂(ℤ_l)` | ★**無い**((G2) に従属) |
| G4 | `mod l` 表現 `Gal(K̄/K) → GL₂(𝔽_l)` | ★**無い**((G1) に従属) |
| G5 | 像が `SL₂` を含むこと | ★**無い**。原文 `Theorem 3.8` の結論そのもの |

★★★**G5 が `[GenEll] Theorem 3.8` の主張である**——
`Lemma 3.1`(`SL₂` の群論、我々は**実装済み**)は
その `mod l` 版を出すための道具である。
-/

namespace ABC3.Interface.GaloisRep

open ABC3.Meta WeierstrassCurve NumberField

/-! ## ★★G3 —— `l` 進表現 -/

/-- ★**`σ ∈ Gal(L/K)` の `E(L)` への作用**——mathlib の `Point.map` そのもの。

★★これは posit ではない(定義できる)。 -/
noncomputable def galPoint {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L]
    [Algebra K L] (W : WeierstrassCurve K) (σ : L ≃ₐ[K] L) :
    (W.baseChange L).toAffine.Point →+ (W.baseChange L).toAffine.Point :=
  WeierstrassCurve.Affine.Point.map (σ.toAlgHom)

/-- ★★**`σ` の `T_l E` への作用**——捩れを保ち、塔と両立するので逆極限に降りる。

★★★これも posit ではない。 -/
noncomputable def galTate {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L]
    [Algebra K L] (W : WeierstrassCurve K) (l : ℕ) (σ : L ≃ₐ[K] L) :
    tateModule (W.baseChange L) l →+ tateModule (W.baseChange L) l where
  toFun := fun f =>
    ⟨fun n => ⟨galPoint W σ ((f : ∀ m : ℕ, torsionPoints (W.baseChange L) (l ^ m)) n), by
      have h : (l ^ n) • (((f : ∀ m : ℕ, torsionPoints (W.baseChange L) (l ^ m)) n :
          (W.baseChange L).toAffine.Point)) = 0 :=
        ((f : ∀ m : ℕ, torsionPoints (W.baseChange L) (l ^ m)) n).2
      show (l ^ n) • (galPoint W σ _) = 0
      rw [← map_nsmul, h, map_zero]⟩, by
    intro n
    show l • (galPoint W σ _) = galPoint W σ _
    rw [← map_nsmul]
    exact congrArg (galPoint W σ) (f.2 n)⟩
  map_zero' := by
    refine Subtype.ext (funext fun n => Subtype.ext ?_)
    exact map_zero _
  map_add' := by
    intro x y
    refine Subtype.ext (funext fun n => Subtype.ext ?_)
    exact map_add _ _ _

/-- **(G3)** `l` 進 Galois 表現 `ρ_{E,l} : Gal(K̄/K) → GL₂(ℤ_l)`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★**原文が名指ししているのはこれである。**

## ★★★★2026-08-20 の訂正(§9-403)

それ以前の版は次の 2 つの欠陥を持っていた。

1. **`rep` が本物の Galois 作用に縛られていなかった**——
   `rep` は素の関数の posit で、条件は `rep_mul`(準同型)と
   `det_rep`(`det_eq_cyclotomic` は同時に posit される)だけだった。
   ★したがって**作用と無関係な準同型で埋まる**。
2. **`det_surjective` は偽であった**——
   `K = ℚ(ζ_l)` を取ると、`σ ∈ Gal(K̄/K)` は `ζ_l` を固定するので
   円分指標の像は `1 + lℤ_l` に含まれ、`ℤ_l^×`(`l ≥ 3` で指数 `l−1`)を
   覆わない。★★数体一般では円分指標は**有限指数の開部分群**にしか落ちない。

★★★訂正:(1) `rep` が `T_l E` への**本物の作用**を表現していることを課し、
(2) 偽の `det_surjective` を落として、代わりに
**`det ρ` が 1 の `l` 冪根に円分指標として作用すること**を課す(= Weil 対の内容、真)。 -/
structure GaloisRepData where
  /-- 台となる Tate 加群。 -/
  toTateModuleData : TateModuleData
  /-- `ρ_{E,l}`。 -/
  rep : {K L : Type} → [Field K] → [DecidableEq K] → [Field L] → [Algebra K L] →
    (W : WeierstrassCurve K) → (l : ℕ) → [Fact l.Prime] → (L ≃ₐ[K] L) →
    Matrix.GeneralLinearGroup (Fin 2) ℤ_[l]
  /-- ★群準同型であること。 -/
  rep_mul : ∀ {K L : Type} [Field K] [DecidableEq K] [Field L] [Algebra K L]
    (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime] (σ τ : L ≃ₐ[K] L),
    rep W l (σ * τ) = rep W l σ * rep W l τ
  /-- ★★★★**`rep` は `T_l E` への本物の作用の行列表示である**。

  ★★これが無いと `rep` は作用と無関係でよい——空虚に埋まる。 -/
  rep_action : ∀ (K L : Type) [Field K] [DecidableEq K] [CharZero K] [Field L] [DecidableEq L]
    [Algebra K L] [IsAlgClosed L] (W : WeierstrassCurve K), W.IsElliptic →
    ∀ (l : ℕ), ∀ _ : Fact l.Prime,
    ∃ e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l]),
      ∀ (σ : L ≃ₐ[K] L) (x : tateModule (W.baseChange L) l),
        e (galTate W l σ x)
          = Matrix.mulVec (rep W l σ : Matrix (Fin 2) (Fin 2) ℤ_[l]) (e x)
  /-- ★★★★**行列式は円分指標である**(Weil 対の内容)。

  `det ρ(σ)` は 1 の `l^n` 乗根への `σ` の作用の指数である。

  ★★これは真であり、かつ **`rep := 1` を殺す**——
  `K = ℚ`、`l ≥ 3` なら `ζ_l` を動かす `σ` があるので `det ρ(σ) ≠ 1`。 -/
  det_cyclotomic : ∀ (K L : Type) [Field K] [DecidableEq K] [CharZero K] [Field L] [DecidableEq L]
    [Algebra K L] [IsAlgClosed L] (W : WeierstrassCurve K), W.IsElliptic →
    ∀ (l : ℕ), ∀ _ : Fact l.Prime, ∀ (σ : L ≃ₐ[K] L) (n : ℕ) (ζ : L), ζ ^ (l ^ n) = 1 →
      σ ζ = ζ ^ ((PadicInt.toZModPow n
        ((rep W l σ : Matrix (Fin 2) (Fin 2) ℤ_[l]).det)).val)

def GaloisRepData.waiting : WaitingFor :=
  { what := "(G3) l 進 Galois 表現 rho_{E,l} : Gal(K̄/K) → GL_2(Z_l) が T_l E への本物の作用の行列表示であることと、その行列式が円分指標であること"
    trackB := "Found/GaloisRep — ★作用そのものは mathlib の `Affine.Point.map`(関手性つき)で**定義済み**(`galTate`)。★★残るのは (i) 作用が Z_l 線型で行列に書けること((G2) の基底を層ごとに使う)と、(ii) **`det = 円分指標`**——これには **Weil 対**が要る(mathlib に 0 件、2026-08-17 実測)" }

/-! ## ★★G4 —— `mod l` 表現 -/

/-- **(G4)** `mod l` 表現 `Gal(K̄/K) → GL₂(𝔽_l)`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★`Lemma 3.1`(`SL₂` の群論)が効くのはこちらである。 -/
structure ModLRepData where
  /-- 台となる `l` 進表現。 -/
  toGaloisRepData : GaloisRepData
  /-- `mod l` 表現。 -/
  repMod : {K L : Type} → [Field K] → [DecidableEq K] → [Field L] → [Algebra K L] →
    (W : WeierstrassCurve K) → (l : ℕ) → [Fact l.Prime] → (L ≃ₐ[K] L) →
    Matrix.GeneralLinearGroup (Fin 2) (ZMod l)
  /-- ★群準同型であること。 -/
  repMod_mul : ∀ {K L : Type} [Field K] [DecidableEq K] [Field L] [Algebra K L]
    (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime] (σ τ : L ≃ₐ[K] L),
    repMod W l (σ * τ) = repMod W l σ * repMod W l τ
  /-- ★★`l` 進表現の還元であること。 -/
  repMod_eq_reduction : ∀ {K L : Type} [Field K] [DecidableEq K] [Field L] [Algebra K L]
    (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime] (σ : L ≃ₐ[K] L),
    (repMod W l σ : Matrix (Fin 2) (Fin 2) (ZMod l))
      = (toGaloisRepData.rep W l σ : Matrix (Fin 2) (Fin 2) ℤ_[l]).map
          (fun a => (PadicInt.toZMod a : ZMod l))

def ModLRepData.waiting : WaitingFor :=
  { what := "(G4) mod l 表現 Gal(K̄/K) → GL_2(F_l) と、それが l 進表現の還元であること"
    trackB := "Found/GaloisRep — ★(G1)(G3) に従属する。★★`PadicInt.toZMod` は mathlib にあるので還元写像は書ける。★★★`Lemma 3.1`(SL_2 の群論)は **我々が Found/GenEll/Lemma31.lean に実装済み**(sorry 0)——効くのはこの mod l の側である" }

/-! ## ★★★G5 —— 像が `SL₂` を含むこと(`Theorem 3.8` の結論) -/

/-- **(G5)** Galois 表現の像が `SL₂(ℤ_l)` を含むこと。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★**これが `[GenEll] Theorem 3.8` の結論そのものである。** -/
structure FullImageData where
  /-- 台となる `mod l` 表現。 -/
  toModLRepData : ModLRepData
  /-- 像が `SL₂` を含むという述語。 -/
  ImageContainsSL2 : {K L : Type} → [Field K] → [DecidableEq K] → [Field L] → [Algebra K L] →
    WeierstrassCurve K → ℕ → Prop
  /-- ★★**述語の内容**——`SL₂(ℤ_l)` のどの元も像に入ること。 -/
  imageContainsSL2_iff : ∀ {K L : Type} [Field K] [DecidableEq K] [Field L] [Algebra K L]
    (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime],
    @ImageContainsSL2 K L _ _ _ _ W l ↔
      ∀ A : Matrix.GeneralLinearGroup (Fin 2) ℤ_[l],
        (A : Matrix (Fin 2) (Fin 2) ℤ_[l]).det = 1 →
        ∃ σ : L ≃ₐ[K] L, toModLRepData.toGaloisRepData.rep W l σ = A

def FullImageData.waiting : WaitingFor :=
  { what := "(G5) Galois 表現の像が SL_2(Z_l) を含むこと —— [GenEll] Theorem 3.8 の結論そのもの"
    trackB := "Found/GaloisRep — ★(G3)(G4) に従属する。★★原文の証明は `Lemma 3.1`(**我々は実装済み**)・`Lemma 3.2`(Tate 曲線、= G6)・`Lemma 3.5` / `Lemma 3.7`(Faltings 高さ、= G8)を使う。★★★したがって S3 の最後の段であり、G1-G4 と G6-G8 がすべて揃ってから着手する" }

/-! ## ★出典の紐付け(`.src`) -/

def GaloisRepData.src : Source :=
  { paper := "GenEll", pdfPage := 19, item := "Theorem 3.8(l 進表現の構成のみ——像の主張は含まない)",
    sectionId := "genell-thm-3-8" }

def ModLRepData.src : Source :=
  { paper := "GenEll", pdfPage := 19, item := "Theorem 3.8(mod l 表現のみ——像の主張は含まない)",
    sectionId := "genell-thm-3-8" }

def FullImageData.src : Source :=
  { paper := "GenEll", pdfPage := 19, item := "Theorem 3.8(像が SL₂ を含むという結論の定式化のみ)",
    sectionId := "genell-thm-3-8" }

end ABC3.Interface.GaloisRep
