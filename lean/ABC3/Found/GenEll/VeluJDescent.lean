/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluDescent
import ABC3.Found.GenEll.PointDescent
import ABC3.Found.GenEll.VeluJExpNeg
import ABC3.Found.GaloisRep.TowerInstances
import ABC3.Found.GaloisRep.DegInfBaseChange
import Mathlib.AlgebraicGeometry.EllipticCurve.IsomOfJ
import Mathlib.AlgebraicGeometry.EllipticCurve.ModelsWithJ
import ABC3.Meta.Claim

/-!
# 第 1444 ブロック —— **同じ `j` の曲線は有限次拡大で同型**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——対偶の道（道 C）の N1

`Skeleton/GenEll/VeluSemistable.lean` に残る `sorry` を対偶で攻める道は

    jExp p E′ < 0
      → E′ の `l`-巡回商は jExp < 0          …… ①
      → その商が E と同じ j を持つ（双対同種） …… ②
      → jExp p E < 0 で矛盾

であり、②は第 1442（`VeluDualJ.lean`）で、①の下位 (i) は第 1443
（`VeluJExpNeg.lean`）で閉じた。★本ファイルは①の残り
——**「`j` が同じなら曲線は同型」を数体へ降ろす段**である。

## ★mathlib に在るもの／無いもの

mathlib は

    WeierstrassCurve.exists_variableChange_of_j_eq (heq : E.j = E'.j) :
        ∃ C : VariableChange F, C • E = E'

を持つが、★**`[IsSepClosed F]` が必須**で、数体 `L` の上ではそのままでは使えない。
☆本ファイルは `AlgebraicClosure L`（`IsSepClosed` のインスタンスは自動で出る）で
`C = ⟨u, r, s, t⟩` を取り、**その 4 成分を含む有限次中間体 `M`** へ降ろす。

★降ろす道具は既に在庫にある:

| 段 | 材料 | 第 |
|---|---|---|
| 有限個の座標は有限次拡大に入る | `exists_finite_subextension` | 1195 |
| 単射な `map` は曲線の等式を落とす | `weierstrassCurve_map_injective` | 1196 |
| `(W ⊗ M) ⊗ L̄ = W ⊗ L̄` | `baseChange_map_intermediate` | 1207 |
| `(C.map φ) • (W.map φ) = (C • W).map φ` | mathlib `map_variableChange` | —— |

☆`u : L̄ˣ` は単元だが、`M` は体なので `Units.mk0` で `Mˣ` に戻せる
——**`u⁻¹` を別に `M` へ入れる必要はない**。

## ★★得られるもの

| 定理 | 内容 |
|---|---|
| `exists_finite_extension_variableChange_of_j_eq` | ★★★★★★**N1 本体**（`j` が同じなら有限次拡大で同型） |
| `exists_finite_extension_variableChange_ofJ` | ★★標準モデル `ofJ E.j` との同型 |
| `exists_numberField_extension_variableChange_of_j_eq` | ★★`NumberField M` つき・`map` の形 |
| `exists_finite_extension_semistableAt_of_jExp_neg` | ★★★★★★**潜在的乗法還元は有限次拡大で乗法還元** |

## ★★★①の「新しい節点」が埋まったこと

`VeluDualJ.lean` の末尾は

> ★①には「潜在的乗法還元は 2 次拡大で乗法還元になる」…… これは新しい節点である。

と書いていた。★`exists_finite_extension_semistableAt_of_jExp_neg` がそれである
——`jExp p E < 0` なら、`L` の**有限次拡大 `M`** と `p` の上のどの素点 `P` でも

    SemistableAt P (E ⊗ M) ∧ jExp P (E ⊗ M) < 0

☆道は 3 本の在庫の合成だけである:

1. `semistableAt_ofJ_j_of_jExp_neg`（第 1443）——`ofJ E.j` は `L` の上で既に半安定
2. `semistableAt_baseChange`（在庫）——半安定は底変換で保たれる
3. 本ファイルの N1 ＋ `semistableAt_variableChange`（在庫）——`M` の上で `E ≅ ofJ E.j`

★★**2 次拡大とは限らない**（`[M : L]` の上界は主張していない）が、
後続が要るのは「有限次」だけである。

## ☆退化していないことの確認

* `SemistableAt` は `Δ = 0` で恒真だが、ここでは `[E.IsElliptic]` を仮定しており
  `Δ` は単元なので `Δ ≠ 0` である。
* `p ∤ l` のような素数の仮定は**一切足していない**（`l` が現れない）。
* 核の座標の整性のような仮定も足していない。
* `jExp` は `j = 0` で `0` を返すが、`jExp p E < 0` は `E.j ≠ 0` を含意し
  （第 1443 の `semistableAt_ofJ_j_of_jExp_neg` が内部で分岐している）、
  本ファイルでは `j` の等式を `ofJ_j` 経由でしか使わないので分岐は不要である。

## ☆逸脱の記録

無し。原典 (GenEll p.17 Lemma 3.5) が「同種なので自動」と畳んだ箇所を埋める
補助節点であり、前提の追加・読み替えはしていない。
★ただし `exists_finite_extension_variableChange_of_j_eq` は
`[NumberField L]` ではなく **`[CharZero L]`** で述べた
（在庫の `exists_finite_subextension` が `[CharZero L]` で足りるため）。
数体はその特別な場合である。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve IsDedekindDomain NumberField ABC3.Meta
open ABC3.Found.GaloisRep

variable {L : Type} [Field L] [CharZero L]

local notation "Lbar" => AlgebraicClosure L

/-! ## ★★★★★★N1 本体——`j` が同じなら有限次拡大で同型 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**同じ `j` の曲線は有限次拡大で同型**——★**無条件**（第 1444）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

`E.j = E'.j` なら、`L` の**有限次拡大 `M ⊂ L̄`** と変数変換
`C : VariableChange M` があって `C • (E ⊗ M) = E' ⊗ M`。

★mathlib の `exists_variableChange_of_j_eq` は `[IsSepClosed F]` を要求するので
数体の上では使えない。☆`L̄` で取った `C = ⟨u, r, s, t⟩` の 4 成分を含む
中間体 `M ≔ L(u, r, s, t)` は `L` 上有限次（第 1195）であり、
`M` は体だから `u` はそこで単元に戻る。
★降ろす一行は第 1196（`map` の単射性）である。 -/
theorem exists_finite_extension_variableChange_of_j_eq
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic] (heq : E.j = E'.j) :
    ∃ M : IntermediateField L Lbar, FiniteDimensional L M ∧
      ∃ C : VariableChange (M : Type),
        C • (E.baseChange (M : Type)) = E'.baseChange (M : Type) := by
  classical
  haveI hEbar : (E.baseChange Lbar).IsElliptic := by
    show (E.map (algebraMap L Lbar)).IsElliptic; infer_instance
  haveI hE'bar : (E'.baseChange Lbar).IsElliptic := by
    show (E'.map (algebraMap L Lbar)).IsElliptic; infer_instance
  have hjbar : (E.baseChange Lbar).j = (E'.baseChange Lbar).j := by
    show (E.map (algebraMap L Lbar)).j = (E'.map (algebraMap L Lbar)).j
    rw [WeierstrassCurve.map_j, WeierstrassCurve.map_j, heq]
  obtain ⟨C, hC⟩ :=
    WeierstrassCurve.exists_variableChange_of_j_eq (E.baseChange Lbar) (E'.baseChange Lbar) hjbar
  obtain ⟨u, r, s, t⟩ := C
  obtain ⟨M, hfin, hmem⟩ := exists_finite_subextension (L := L)
    ({((u : Lbar), r), (s, t)} : Finset (Lbar × Lbar))
  obtain ⟨hu, hr⟩ := hmem ((u : Lbar), r) (by simp)
  obtain ⟨hs, ht⟩ := hmem (s, t) (by simp)
  have hu0 : (⟨(u : Lbar), hu⟩ : (M : Type)) ≠ 0 := fun h => u.ne_zero (congrArg Subtype.val h)
  refine ⟨M, hfin, ⟨Units.mk0 ⟨(u : Lbar), hu⟩ hu0, ⟨r, hr⟩, ⟨s, hs⟩, ⟨t, ht⟩⟩, ?_⟩
  refine weierstrassCurve_map_injective (algebraMap (M : Type) Lbar)
    (algebraMap (M : Type) Lbar).injective ?_
  rw [← WeierstrassCurve.map_variableChange]
  have hCM : (⟨Units.mk0 ⟨(u : Lbar), hu⟩ hu0, ⟨r, hr⟩, ⟨s, hs⟩, ⟨t, ht⟩⟩ :
      VariableChange (M : Type)).map (algebraMap (M : Type) Lbar)
      = (⟨u, r, s, t⟩ : VariableChange Lbar) := by
    ext <;> rfl
  rw [hCM, baseChange_map_intermediate, baseChange_map_intermediate]
  exact hC

/-- ★★**標準モデルとの同型**——`E` は有限次拡大の上で `ofJ E.j` と同型である。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`ofJ_j : (ofJ j).j = j` を N1 に入れるだけである。 -/
theorem exists_finite_extension_variableChange_ofJ [DecidableEq L]
    (E : WeierstrassCurve L) [E.IsElliptic] :
    ∃ M : IntermediateField L Lbar, FiniteDimensional L M ∧
      ∃ C : VariableChange (M : Type),
        C • (E.baseChange (M : Type)) = (WeierstrassCurve.ofJ E.j).baseChange (M : Type) :=
  exists_finite_extension_variableChange_of_j_eq E (WeierstrassCurve.ofJ E.j)
    (WeierstrassCurve.ofJ_j E.j).symm

/-- ★★**N1 の数体版**——`M` が数体であることと、`baseChange` を `map` で書いた形。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`NumberField M` は第 1222（`numberField_and_towers`）から出る。 -/
theorem exists_numberField_extension_variableChange_of_j_eq
    (K : Type) [Field K] [NumberField K] (E E' : WeierstrassCurve K)
    [E.IsElliptic] [E'.IsElliptic] (heq : E.j = E'.j) :
    ∃ M : IntermediateField K (AlgebraicClosure K), FiniteDimensional K M ∧
      NumberField (M : Type) ∧
      ∃ C : VariableChange (M : Type),
        C • (E.map (algebraMap K (M : Type))) = E'.map (algebraMap K (M : Type)) := by
  obtain ⟨M, hfin, C, hC⟩ := exists_finite_extension_variableChange_of_j_eq E E' heq
  exact ⟨M, hfin, (numberField_and_towers K M).1, C, hC⟩

/-! ## ★★★★★★潜在的乗法還元は有限次拡大で乗法還元 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★
**潜在的乗法還元は有限次拡大で乗法還元になる**——★**無条件**（第 1444）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

`jExp p E < 0`（＝ `v_p(j) < 0`）なら、`L` の**有限次拡大 `M`** があって、
`p` の上にあるどの素点 `P` でも

    SemistableAt P (E ⊗ M)   かつ   jExp P (E ⊗ M) < 0

★★これが `VeluDualJ.lean` が「新しい節点である」と書いた段である。
☆`E` そのものは加法還元でありうる（2 次のねじれ）ので体の拡大が要るが、
`j` だけで決まる標準モデル `ofJ E.j` は第 1443 で**`L` の上で既に**半安定だから、
必要なのは「`E` と `ofJ E.j` が `M` の上で同型」——すなわち N1 である。

★`jExp` の側は `jExp_baseChange`（`v_P(j) = e·v_p(j)`）と `e ≠ 0` から出る。 -/
theorem exists_finite_extension_semistableAt_of_jExp_neg
    (K : Type) [Field K] [NumberField K] [DecidableEq K]
    (p : HeightOneSpectrum (𝓞 K)) (E : WeierstrassCurve K) [E.IsElliptic]
    (hj : jExp p E < 0) :
    ∃ (M : IntermediateField K (AlgebraicClosure K)) (hfin : FiniteDimensional K M),
      haveI : FiniteDimensional K (M : Type) := hfin
      letI : NumberField (M : Type) := NumberField.of_module_finite K M
      letI : IsScalarTower (𝓞 K) K (M : Type) := isScalarTower_ringOfIntegers_base K M
      letI : IsScalarTower (𝓞 K) (𝓞 (M : Type)) (M : Type) := isScalarTower_ringOfIntegers_top K M
      ∀ P : HeightOneSpectrum (𝓞 (M : Type)), P.asIdeal.LiesOver p.asIdeal →
        SemistableAt P (E.baseChange (M : Type)) ∧ jExp P (E.baseChange (M : Type)) < 0 := by
  obtain ⟨M, hfin, C, hC⟩ := exists_finite_extension_variableChange_ofJ (L := K) E
  refine ⟨M, hfin, ?_⟩
  haveI : FiniteDimensional K (M : Type) := hfin
  letI : NumberField (M : Type) := NumberField.of_module_finite K M
  letI : IsScalarTower (𝓞 K) K (M : Type) := isScalarTower_ringOfIntegers_base K M
  letI : IsScalarTower (𝓞 K) (𝓞 (M : Type)) (M : Type) := isScalarTower_ringOfIntegers_top K M
  intro P hP
  haveI := hP
  haveI : ((WeierstrassCurve.ofJ E.j).baseChange (M : Type)).IsElliptic := by
    show ((WeierstrassCurve.ofJ E.j).map (algebraMap K (M : Type))).IsElliptic
    infer_instance
  have hss : SemistableAt p (WeierstrassCurve.ofJ E.j) :=
    (semistableAt_ofJ_j_of_jExp_neg p E hj).1
  have hssM : SemistableAt P ((WeierstrassCurve.ofJ E.j).baseChange (M : Type)) :=
    semistableAt_baseChange K (M : Type) p P (WeierstrassCurve.ofJ E.j) hss
  have hEq : E.baseChange (M : Type)
      = C⁻¹ • ((WeierstrassCurve.ofJ E.j).baseChange (M : Type)) := by
    rw [← hC, inv_smul_smul]
  refine ⟨?_, ?_⟩
  · rw [hEq]
    exact semistableAt_variableChange P _ _ hssM
  · rw [jExp_baseChange K (M : Type) p P E]
    exact mul_neg_of_pos_of_neg
      (by exact_mod_cast Nat.pos_of_ne_zero (ramificationIdx_ne_zero p P)) hj

/-! ## ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def exists_finite_extension_variableChange_of_j_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(同じ j の曲線は有限次拡大で同型。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_finite_extension_variableChange_of_j_eq.needs : List ProofObligation :=
  [ .citation "[mathlib]" "exists_variableChange_of_j_eq(分離閉体の上での同型)"
      (.inMathlib "WeierstrassCurve.exists_variableChange_of_j_eq") 1,
    .citation "[ABC3]" "exists_finite_subextension(第 1195、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_finite_subextension") 1,
    .citation "[ABC3]" "weierstrassCurve_map_injective(第 1196、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.weierstrassCurve_map_injective") 1,
    .implicitStep
      ("★★★★★**2026-09-06（第 1444）**——mathlib の `exists_variableChange_of_j_eq` は " ++
       "`[IsSepClosed F]` を要求するので数体では使えない。" ++
       "☆`L̄`（`IsSepClosed` は自動）で変数変換 `C = ⟨u, r, s, t⟩` を取り、" ++
       "その 4 成分を含む中間体 `L(u, r, s, t)` へ降ろす。" ++
       "★`u` は `M` が体だから `Units.mk0` で単元に戻り、`u⁻¹` を別に足す必要はない。") 17 ]

def exists_finite_extension_variableChange_ofJ.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(有限次拡大の上で標準モデル ofJ E.j と同型。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_numberField_extension_variableChange_of_j_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(同じ j の曲線は有限次拡大で同型——数体版。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_finite_extension_semistableAt_of_jExp_neg.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(潜在的乗法還元は有限次拡大で乗法還元になる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_finite_extension_semistableAt_of_jExp_neg.needs : List ProofObligation :=
  [ .citation "[ABC3]" "semistableAt_ofJ_j_of_jExp_neg(第 1443、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.semistableAt_ofJ_j_of_jExp_neg") 1,
    .citation "[ABC3]" "semistableAt_baseChange(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.semistableAt_baseChange") 1,
    .citation "[ABC3]" "semistableAt_variableChange(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.semistableAt_variableChange") 1,
    .implicitStep
      ("★★★★★★**2026-09-06（第 1444）**——`VeluDualJ.lean` が" ++
       "「①には潜在的乗法還元は有限次拡大で乗法還元になることが要る。" ++
       "これは新しい節点である」と書いた段である。" ++
       "☆`E` 自身は加法還元でありうるが、`j` だけで決まる標準モデル `ofJ E.j` は" ++
       "第 1443 で `L` の上で既に半安定なので、要るのは" ++
       "「`E` と `ofJ E.j` が有限次拡大の上で同型」——すなわち本ファイルの N1 だけである。") 17 ]

end ABC3.Found.GenEll
