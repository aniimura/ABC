/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.DivisorDescent
import ABC3.Found.GenEll.CartierPullback
import ABC3.Found.GenEll.BDClass

/-!
# [GenEll] Remark 1.5.1 —— **`Σ` の外での導手の一致**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

## ★★★★★★★★★★これで `Remark 1.5.1` の 5 段が揃う

| 段 | 内容 | どこ |
|---|---|---|
| 1 | 生成ファイバーの同型が `ℤ[Σ⁻¹]` へ延びる | `IsoDescent.lean` |
| 2 | 因子 `D` も一緒に延びる | `DivisorDescent.lean` |
| 3 | ★**`Σ` の外で導手が一致する** | ★**本ファイル** |
| 4 | `Σ` の上の寄与の一様上界 | `LogCondSigma.lean` |
| 5 | `log-diff` は `X_ℚ` だけに依る | `LogDiff.lean` |

## ★★★★★★機構は 2 つに割れる

**算術の側**（`conductorADiv_eq_of_comap_eq`）——
**点を動かすか因子を動かすかは同じ**である:

    pullbackIdeal F D (x ≫ p) = pullbackIdeal F (D.comap p) x

★これは `Scheme.IdealSheafData.comap_comp` の言い換えにすぎない。
★★したがって「因子が段 `n` で一致する」ことさえ言えれば、導手は**点ごとに**一致する。

**幾何の側**（`comap_eq_of_square`）——
`DivisorDescent.lean` が出す可換正方形を**引き戻し正方形**に読み替える:

* `bcBC` は `iZ` の底変換だから、`bcBC ≫ p = snd ≫ iZ` は**引き戻し**である
  （`bcBC_isPullback`、`IsPullback.of_right` で貼る）。
* 両側が同型の可換正方形も**引き戻し**である（`IsPullback.of_horiz_isIso`）。
* 引き戻し正方形の核は `comap` である（`ker_eq_comap_of_isPullback`）。

★★★**3 つを繋ぐと `iZ.ker.comap p = (iZ′.ker.comap p′).comap φ` が出る。**

## ★★★★ここでも mathlib の欠落を迂回している

イデアル層の水準で直接やると
`(I.comap f).ideal U = (I.ideal V).map (f.appLE V U)`（mathlib に無い）が要る。
★**核と引き戻し正方形だけで書けば要らない**——`SectionDescent.lean` の測定を見よ。
-/

namespace ABC3.Found.GenEll

open CategoryTheory Limits AlgebraicGeometry

/-! ## ★★★★★★算術の側 —— 点を動かすか因子を動かすか -/

variable (F : Type) [Field F] [NumberField F]

omit [NumberField F] in
/-- ★★★★**点を動かすか因子を動かすかは同じ**。

原文 (GenEll p.8):
> Then the morphism Spec(OF ) → X determined by x [where we recall that X is proper!] allows one to pull-back the divisor D to Spec(OF ) so as to obtain an eﬀective divisor Dx on Spec(OF ).

★中身は `Scheme.IdealSheafData.comap_comp` の言い換えである。 -/
theorem pullbackIdeal_target {X X' : Scheme} (D' : X'.IdealSheafData) (φ : X ⟶ X')
    (xF : specRingOfIntegers F ⟶ X) :
    pullbackIdeal F D' (xF ≫ φ) = pullbackIdeal F (D'.comap φ) xF := by
  simp [pullbackIdeal, Scheme.IdealSheafData.comap_comp]

/-- ★★★★★★**因子が段 `n` で一致していれば、導手は点ごとに一致する**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

★仮定 `h` は「`D` と `D′` が `W` の上で（`φ` を通して）同じ因子になる」こと。 -/
theorem conductorADiv_eq_of_comap_eq {W W' X X' : Scheme}
    (D : X.IdealSheafData) (D' : X'.IdealSheafData)
    (p : W ⟶ X) (p' : W' ⟶ X') (φ : W ⟶ W')
    (h : D.comap p = (D'.comap p').comap φ)
    (xF : specRingOfIntegers F ⟶ W) :
    conductorADiv F D (xF ≫ p) = conductorADiv F D' (xF ≫ φ ≫ p') := by
  unfold conductorADiv
  rw [pullbackIdeal_target, pullbackIdeal_target, h, Scheme.IdealSheafData.comap_comp]

/-- ★★★★★★同じことを `log-cond` の水準で。 -/
theorem logCond_eq_of_comap_eq {W W' X X' : Scheme}
    (D : X.IdealSheafData) (D' : X'.IdealSheafData)
    (p : W ⟶ X) (p' : W' ⟶ X') (φ : W ⟶ W')
    (h : D.comap p = (D'.comap p').comap φ)
    (xF : specRingOfIntegers F ⟶ W) :
    logCond F D (xF ≫ p) = logCond F D' (xF ≫ φ ≫ p') := by
  unfold logCond
  rw [conductorADiv_eq_of_comap_eq F D D' p p' φ h xF]

/-! ## ★★★★★★★★幾何の側 —— 引き戻し正方形の核 -/

/-- ★★★★**引き戻し正方形なら、核は `comap` である**。

★`mathlib` の `ker_fst_of_isClosedImmersion` は**標準の引き戻し**についてのもの。
一般の引き戻し正方形へは `IsPullback.isoPullback` で移す（同型に沿った核は不変）。 -/
theorem ker_eq_comap_of_isPullback {P Z W X : Scheme.{0}}
    (u : P ⟶ W) (v : P ⟶ Z) (p : W ⟶ X) (iZ : Z ⟶ X)
    [IsClosedImmersion iZ] (H : IsPullback u v p iZ) :
    u.ker = iZ.ker.comap p := by
  rw [← Scheme.IdealSheafData.ker_fst_of_isClosedImmersion iZ p, ← H.isoPullback_hom_fst,
    Scheme.Hom.ker_comp_of_isIso]

/-- ★★★★★★**底変換 `bcBC` は引き戻し正方形をなす**。

    X_n ─bcBC→ ... と Z ─iZ→ X の正方形

★`Spec ℤ` は終対象なので `Z` の構造射は `iZ ≫ f` に他ならない。
★★機構は横方向の貼り合わせ（`IsPullback.of_right`）で、
右の正方形が `bcObj f n` の定義そのもの、外側が `bcObj (iZ ≫ f) n` の定義そのものである。 -/
theorem bcBC_isPullback {Z X : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (iZ : Z ⟶ X) (n : ℕᵒᵖ) :
    IsPullback (bcBC (iZ ≫ f) f iZ n)
      (pullback.snd (overRatTowerDiagram.obj n).hom (iZ ≫ f))
      (pullback.snd (overRatTowerDiagram.obj n).hom f) iZ := by
  refine IsPullback.of_right ?_ (bcBC_snd (iZ ≫ f) f iZ n)
    (IsPullback.of_hasPullback (overRatTowerDiagram.obj n).hom f)
  rw [bcBC_fst]
  exact IsPullback.of_hasPullback (overRatTowerDiagram.obj n).hom (iZ ≫ f)

/-- ★★閉埋め込みの底変換は閉埋め込み。 -/
instance bcBC_isClosedImmersion {Z X : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (iZ : Z ⟶ X) [IsClosedImmersion iZ] (n : ℕᵒᵖ) :
    IsClosedImmersion (bcBC (iZ ≫ f) f iZ n) :=
  MorphismProperty.of_isPullback (bcBC_isPullback f iZ n).flip ‹_›

/-- ★★★★★★★★**可換正方形＋両側同型 ⟹ 因子が段 `n` で一致する**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

★これが `DivisorDescent.lean` の `exists_pair_descent` が出す正方形の**読み替え**である。
★★機構は 3 つ: (a) `bcBC` の正方形は引き戻し、(b) 両側同型の可換正方形も引き戻し、
(c) 引き戻し正方形の核は `comap`。 -/
theorem comap_eq_of_square {Z Z' X X' : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ))
    (iZ : Z ⟶ X) (iZ' : Z' ⟶ X') [IsClosedImmersion iZ] [IsClosedImmersion iZ']
    {n : ℕᵒᵖ}
    (φ : bcObj f n ⟶ bcObj f' n) (ψ : bcObj (iZ ≫ f) n ⟶ bcObj (iZ' ≫ f') n)
    [IsIso φ] [IsIso ψ]
    (hsq : ψ ≫ bcBC (iZ' ≫ f') f' iZ' n = bcBC (iZ ≫ f) f iZ n ≫ φ) :
    iZ.ker.comap (pullback.snd (overRatTowerDiagram.obj n).hom f)
      = (iZ'.ker.comap (pullback.snd (overRatTowerDiagram.obj n).hom f')).comap φ := by
  have h1 := ker_eq_comap_of_isPullback _ _ _ _ (bcBC_isPullback f iZ n)
  have h2 := ker_eq_comap_of_isPullback _ _ _ _ (bcBC_isPullback f' iZ' n)
  have hpb : IsPullback (bcBC (iZ ≫ f) f iZ n) ψ φ (bcBC (iZ' ≫ f') f' iZ' n) :=
    (IsPullback.of_horiz_isIso ⟨hsq⟩).flip
  have h3 := ker_eq_comap_of_isPullback _ _ _ _ hpb
  rw [← h1, h3, ← h2]

/-! ## ★★★★★★★★★★到達点 —— 導手の一致 -/

/-- ★★★★★★★★★★**[GenEll] Remark 1.5.1 —— 段 3: 導手の一致**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

`DivisorDescent.lean` の `exists_pair_descent` が出す
「有限段 `n` での同型＋可換正方形」から、
**`X ×_ℤ ℤ[1/n!]` の点ごとに導手が一致する**ことが出る。

★★**「`Σ` の外で」という原文の条件は、点が段 `n` を経由することに相当する**
——`n!` が可逆な素点でだけ `X_n` の点になれる。 -/
theorem conductorADiv_eq_of_square {Z Z' X X' : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ))
    (iZ : Z ⟶ X) (iZ' : Z' ⟶ X') [IsClosedImmersion iZ] [IsClosedImmersion iZ']
    {n : ℕᵒᵖ}
    (φ : bcObj f n ⟶ bcObj f' n) (ψ : bcObj (iZ ≫ f) n ⟶ bcObj (iZ' ≫ f') n)
    [IsIso φ] [IsIso ψ]
    (hsq : ψ ≫ bcBC (iZ' ≫ f') f' iZ' n = bcBC (iZ ≫ f) f iZ n ≫ φ)
    (xF : specRingOfIntegers F ⟶ bcObj f n) :
    conductorADiv F iZ.ker (xF ≫ pullback.snd (overRatTowerDiagram.obj n).hom f)
      = conductorADiv F iZ'.ker (xF ≫ φ ≫ pullback.snd (overRatTowerDiagram.obj n).hom f') :=
  conductorADiv_eq_of_comap_eq F _ _ _ _ φ (comap_eq_of_square f f' iZ iZ' φ ψ hsq) xF

/-- ★★★★★★★★★★同じことを `log-cond` の水準で。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

★★**定数差すら出ない**（`C = 0`）——段 `n` の上では因子そのものが一致するからである。 -/
theorem logCond_eq_of_square {Z Z' X X' : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ))
    (iZ : Z ⟶ X) (iZ' : Z' ⟶ X') [IsClosedImmersion iZ] [IsClosedImmersion iZ']
    {n : ℕᵒᵖ}
    (φ : bcObj f n ⟶ bcObj f' n) (ψ : bcObj (iZ ≫ f) n ⟶ bcObj (iZ' ≫ f') n)
    [IsIso φ] [IsIso ψ]
    (hsq : ψ ≫ bcBC (iZ' ≫ f') f' iZ' n = bcBC (iZ ≫ f) f iZ n ≫ φ)
    (xF : specRingOfIntegers F ⟶ bcObj f n) :
    logCond F iZ.ker (xF ≫ pullback.snd (overRatTowerDiagram.obj n).hom f)
      = logCond F iZ'.ker (xF ≫ φ ≫ pullback.snd (overRatTowerDiagram.obj n).hom f') :=
  logCond_eq_of_comap_eq F _ _ _ _ φ (comap_eq_of_square f f' iZ iZ' φ ψ hsq) xF

/-- ★★★★★★★★★★**[GenEll] Remark 1.5.1 —— 段 `n` の点の上では BD-同値**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

★★**定数差すら出ない**（`C = 0`）——段 `n` の上では因子そのものが一致するからである。
★原文が `≈`（BD-同値）で済ませているのは、`X` の点全体を見ると
`Σ` の上の寄与が残るからであり（段 4、`LogCondSigma.lean`）、
**段 `n` の点に限ればその寄与も消える**。 -/
theorem remark_1_5_1_bdeq_of_square {Z Z' X X' : Scheme.{0}}
    (f : X ⟶ Spec (CommRingCat.of ℤ)) (f' : X' ⟶ Spec (CommRingCat.of ℤ))
    (iZ : Z ⟶ X) (iZ' : Z' ⟶ X') [IsClosedImmersion iZ] [IsClosedImmersion iZ']
    {n : ℕᵒᵖ}
    (φ : bcObj f n ⟶ bcObj f' n) (ψ : bcObj (iZ ≫ f) n ⟶ bcObj (iZ' ≫ f') n)
    [IsIso φ] [IsIso ψ]
    (hsq : ψ ≫ bcBC (iZ' ≫ f') f' iZ' n = bcBC (iZ ≫ f) f iZ n ≫ φ) :
    BDeq (fun xF : specRingOfIntegers F ⟶ bcObj f n =>
            logCond F iZ.ker (xF ≫ pullback.snd (overRatTowerDiagram.obj n).hom f))
         (fun xF => logCond F iZ'.ker
            (xF ≫ φ ≫ pullback.snd (overRatTowerDiagram.obj n).hom f')) :=
  ⟨0, fun xF => by
    dsimp only
    rw [logCond_eq_of_square F f f' iZ iZ' φ ψ hsq xF, sub_self, abs_zero]⟩

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` はまだ置かない。** 残っているのは
「`X` の点のうち `Σ` の外の素点でだけ `X_n` の点になれる」という**橋**である
——本ファイルは `X_n` の点について述べており、原文は `X` の点について述べている。 -/

def remark_1_5_1_bdeq_of_square.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(段 n の点の上で log-cond が BD-同値。X の点への橋は含まない)",
    sectionId := "genell-rem-1-5-1" }

def pullbackIdeal_target.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iv)(点を動かすか因子を動かすかは同じ)",
    sectionId := "genell-def-1-5" }

def conductorADiv_eq_of_square.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(段 3——段 n の点について導手が一致する。X の点への橋は含まない)",
    sectionId := "genell-rem-1-5-1" }

def conductorADiv_eq_of_square.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_pair_descent(対 (X, D) の spreading out——正方形を出す)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_pair_descent") 9,
    .citation "[mathlib]" "Scheme.IdealSheafData.ker_fst_of_isClosedImmersion(標準の引き戻しの核)"
      (.inMathlib "AlgebraicGeometry.Scheme.IdealSheafData.ker_fst_of_isClosedImmersion") 9,
    .citation "[mathlib]" "IsPullback.of_right(横方向の貼り合わせ)"
      (.inMathlib "CategoryTheory.IsPullback.of_right") 9,
    .citation "[mathlib]" "IsPullback.of_horiz_isIso(両側同型の可換正方形は引き戻し)"
      (.inMathlib "CategoryTheory.IsPullback.of_horiz_isIso") 9,
    .citation "[mathlib]" "Scheme.IdealSheafData.comap_comp(点を動かすか因子を動かすか)"
      (.inMathlib "AlgebraicGeometry.Scheme.IdealSheafData.comap_comp") 9,
    .implicitStep
      ("★★★★イデアル層の水準で直接やると mathlib に無い ideal_comap が要る。" ++
       "★核と引き戻し正方形だけで書けば要らない" ++
       "(ResearchPaper/mathlib-gap.json の ideal-comap-on-affine-opens)") 9,
    .implicitStep
      ("★★★★★残る橋: 原文は X の点と『Σ の外の素点 v』について述べ、" ++
       "本ファイルは X_n = X ×_ℤ ℤ[1/n!] の点について述べている。" ++
       "★両者を繋ぐのは「v が n! を割らないなら x は局所的に X_n を経由する」" ++
       "という段であり、これは可換環論(局所化)の作業である") 9 ]

end ABC3.Found.GenEll
