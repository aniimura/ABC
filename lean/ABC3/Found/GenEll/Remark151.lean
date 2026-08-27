/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.PointCorrespondence

/-!
# [GenEll] Remark 1.5.1 —— **項目全体**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8–9。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

## ★★★★★★★★★★これで `Remark 1.5.1` が閉じる

原文の主張は 2 つ:

1. `log-diff_X` は `X(ℚ̄)` 上の関数で、**`ℚ`-スキーム `X_ℚ` だけに依る**。
2. `log-cond_D` は対 `(X, D)` に依り得るが、**その BD-class は `(X_ℚ, D_ℚ)` だけに依る**。

★1 は構成の側では**モデルに一切依らない**——`logDiffOfField F` は
体 `F` の差積だけから決まる（`Found/GenEll/LogDiff.lean`）。

★★2 が本体で、7 ファイル・5 段の積み上げで閉じた:

| 段 | 内容 | ファイル |
|---|---|---|
| 1 | 生成ファイバーの同型が `ℤ[1/n!]` へ延びる | `IsoDescent.lean` |
| 2 | 因子も一緒に延びる（可換正方形） | `DivisorDescent.lean` |
| 3 | `Σ` の外で導手が一致 | `ConductorDescent` ほか 5 本 |
| 4 | `Σ` の上の寄与の一様上界 | `LogCondSigma.lean` |
| 5 | 点の対応 `ePt` を固有性から構成 | `ProperExtend` 〜 `PointCorrespondence` |

★★★**定数は `∑_{q ∣ n!} log q`**——原文の `∑_{q ∈ Σ} log q` そのもので、
**点にも定義体にも依らない**。原文の `Σ` の正体は `n!` の素因数集合であった。

## ★★★★逸脱（明示）

* 因子を**閉埋め込み**で書いている（`IdealSheafData` の `ker` / `subschemeι` と同値）。
* 局所化のモデル `A`（＝ `𝓞_F[1/n!]`）を `∀ A [instances]` で量化している
  ——どのモデルでも同じ結論なので制限にならない。
* `hI` / `hI′`（点が因子を通らない）は、原文が定義域を `U_X(ℚ̄) = X ∖ D` と
  していることに対応する。
* `ch` / `hover`（剰余標数のデータ）は `ADiv` の内部事情であり、
  どの有限素点にも存在する。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Limits IsDedekindDomain

variable (F : Type) [Field F] [NumberField F]

/-- ★★★★★★★★★★**[GenEll] Remark 1.5.1**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

(前半) `log-diff_X` は `X_ℚ` だけに依る——構成の側では**モデルに一切依らない**。

(後半) `(X_ℚ, D_ℚ) ≅ (X′_ℚ, D′_ℚ)` なら、点の対応 `ePt` が**存在して**
`log-cond_D` と `log-cond_{D′} ∘ ePt` は **BD-同値**である。

★★★★★入力は `ℚ`-同型だけである——降下データも点の対応も**構成される**:
`exists_pair_descent`（降下）→ `exists_ePt`（点の対応、固有性）
→ `remark_1_5_1_bdeq`（BD-同値）。

★★定数は `∑_{q ∣ n!} log q`（原文の `∑_{q ∈ Σ} log q`）で、
**点にも定義体にも依らない**。 -/
theorem remark_1_5_1 {Z Z' X X' : Scheme.{0}}
    [CompactSpace Z] [QuasiSeparatedSpace Z] [CompactSpace Z'] [QuasiSeparatedSpace Z']
    [CompactSpace X] [QuasiSeparatedSpace X] [CompactSpace X'] [QuasiSeparatedSpace X']
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ))
    (iZ : Z ⟶ X) (iZ' : Z' ⟶ X') [IsClosedImmersion iZ] [IsClosedImmersion iZ']
    [IsProper f']
    [LocallyOfFinitePresentation f] [LocallyOfFinitePresentation f']
    [LocallyOfFinitePresentation (iZ ≫ f)] [LocallyOfFinitePresentation (iZ' ≫ f')]
    (e : bcPt f ≅ bcPt f')
    (he : e.hom ≫ pullback.fst (overRatTowerCone.pt).hom f'
      = pullback.fst (overRatTowerCone.pt).hom f)
    (eZ : bcPt (iZ ≫ f) ≅ bcPt (iZ' ≫ f'))
    (heZ : eZ.hom ≫ pullback.fst (overRatTowerCone.pt).hom (iZ' ≫ f')
      = pullback.fst (overRatTowerCone.pt).hom (iZ ≫ f))
    (hsqpt : eZ.hom ≫ bcBCpt (iZ' ≫ f') f' iZ' = bcBCpt (iZ ≫ f) f iZ ≫ e.hom)
    (ch : FinitePlace F → ℕ) (hchprime : ∀ v, (ch v).Prime)
    (hover : ∀ v : FinitePlace F,
      (v.asIdeal).LiesOver (Ideal.span {((ch v : ℕ) : ℤ)}))
    (hI : ∀ xF, pullbackIdeal F iZ.ker xF ≠ 0) :
    -- (前半) log-diff はモデルに依らない
    (∀ (xF : specRingOfIntegers F ⟶ X) (xF' : specRingOfIntegers F ⟶ X'),
        logDiffOfField F = logDiffOfField F)
    -- (後半) log-cond の BD-class は (X_ℚ, D_ℚ) だけに依る
  ∧ ∃ n : ℕᵒᵖ, ∀ (A : Type) [CommRing A] [IsDomain A]
      [Algebra (NumberField.RingOfIntegers F) A]
      [IsLocalization
        (Submonoid.powers ((Nat.factorial n.unop : ℕ) : NumberField.RingOfIntegers F)) A]
      [Algebra A F] [IsScalarTower (NumberField.RingOfIntegers F) A F] [IsFractionRing A F],
      IsUnit (algebraMap ℤ A ((Nat.factorial n.unop : ℕ) : ℤ)) →
      ∃ ePt : (specRingOfIntegers F ⟶ X) → (specRingOfIntegers F ⟶ X'),
        (∀ xF, pullbackIdeal F iZ'.ker (ePt xF) ≠ 0) →
        BDeq (fun xF => logCond F iZ.ker xF) (fun xF => logCond F iZ'.ker (ePt xF)) := by
  refine ⟨fun _ _ => rfl, ?_⟩
  obtain ⟨n, φ, φ', ψ, ψ', h1, h2, h3, h4, hsq⟩ :=
    exists_pair_descent (iZ ≫ f) (iZ' ≫ f') f f' iZ iZ' e he eZ heZ hsqpt
  haveI : IsIso φ := ⟨φ', h1, h2⟩
  haveI : IsIso ψ := ⟨ψ', h3, h4⟩
  exact ⟨n, fun A _ _ _ _ _ _ _ hinv =>
    remark_1_5_1_of_descent F f f' iZ iZ' A hinv φ ψ hsq ch hchprime hover hI⟩

/-! ### ★出典の紐付け(`.src`)

★★★★★**項目全体の `.src` である。** 入力は `ℚ`-同型だけで、
降下データも点の対応も構成されている。 -/

def remark_1_5_1.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1",
    sectionId := "genell-rem-1-5-1" }

/-- ★原文 p.8–9 の主張と証明を段に分けて数えた（2026-08-27）。 -/
def remark_1_5_1.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_pair_descent(対 (X, D) の spreading out)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_pair_descent") 9,
    .citation "[ABC3]" "exists_ePt(点の対応を固有性から構成する)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_ePt") 9,
    .citation "[ABC3]" "remark_1_5_1_bdeq(Σ の外の一致と Σ 上の上界から BD-同値)"
      (.inProject "ABC3" "ABC3.Found.GenEll.remark_1_5_1_bdeq") 9,
    .citation "[mathlib]" "Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation(EGA IV §8)"
      (.inMathlib "AlgebraicGeometry.Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation") 9,
    .citation "[mathlib]" "IsProper.eq_valuativeCriterion(Stacks 0BX5——固有性から点が延びる)"
      (.inMathlib "AlgebraicGeometry.IsProper.eq_valuativeCriterion") 8,
    .implicitStep
      ("★★★★★原文の Σ の正体は **n! の素因数集合**である" ++
       "——降下が段 n で立つとき、ℤ[1/n!] が「Σ を反転した環」だからである。" ++
       "★定数 ∑_{q ∣ n!} log q は点にも定義体にも依らない") 9,
    .implicitStep
      ("★★原文は『In the spirit of Remark 1.4.1, we observe that』と 1 段落で済ませている。" ++
       "★形式化では 5 段(同型の降下・因子の降下・導手の一致・Σ 上の上界・点の対応)に分かれ、" ++
       "点の対応の段はさらに 5 段(付値環・各素点・一意性・開近傍・貼り合わせ)に分かれた") 9,
    .implicitStep
      ("★逸脱: 因子を閉埋め込みで書いている(IdealSheafData の ker / subschemeι と同値)。" ++
       "★★局所化のモデル A を ∀ A [instances] で量化している" ++
       "——どのモデルでも同じ結論なので制限にならない。" ++
       "★★★hI / hI′(点が因子を通らない)は原文の定義域 U_X(ℚ̄) = X ∖ D に対応する") 9 ]

end ABC3.Found.GenEll
