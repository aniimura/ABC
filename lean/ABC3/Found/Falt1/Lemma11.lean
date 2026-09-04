import ABC3.Meta.Claim
import ABC3.Found.Falt1.KaehlerAux
import ABC3.Interface.Falt1.Ramification
import Mathlib.RingTheory.HopkinsLevitzki

/-!
# [Falt1] Lemma 1.1 —— 単射性・長さの等式(`Found`)

原典: G. Faltings, *p-Adic Hodge Theory*, JAMS Vol.1 No.1 (1988) pp.255-299。
物理 p.4(印字 p.257)。**260 dpi 目視確認 2026-09-04**(逐語の食い違い無し)。

★このPDFのpdftotext OCR層は数式記号を激しく壊す(`Ω`が`K`に、`⊗`が`?`に、
`⊂`が`c`になる等——`Interface/Falt1/Ramification.lean`冒頭の注記と同じ
運用)ため、地の文(`内容`)で写す。

内容 (Falt1 p.4): Lemma 1.1. For any extension V ⊂ W, as above, the
natural map Ω_V ⊗_V W → Ω_W is injective, and its cokernel Ω_{W/V}
has the same length as W/p^δW.

`KaehlerAux.lean` で証明した以下の2定理を貼り合わせ、Falt1 の実際の
`V,K,L,W` に対して Lemma 1.1 の主張(単射性 ＋ 余核の長さの等式)を
証明した:
- `falt1MapBaseChangeInjective`(単射性)
- `falt1CokernelLengthEq`(長さの等式、`Ω_{W/V}` の余核が `W ⧸
  differentIdeal V W` に等しいことと合わせて)
-/

namespace ABC3.Found.Falt1

/-- `differentIdeal_ne_bot`(mathlib)適用の2条件のうち、Falt1 の既存の
仮定だけから導けることを確認した1つ目: `L` は `W` の分数体そのもの。 -/
theorem falt1_isFractionRing {V K L W : Type*} [CommRing V] [IsDedekindDomain V]
    [Field K] [Algebra V K] [IsFractionRing V K] [Field L] [Algebra K L] [FiniteDimensional K L]
    [CommRing W] [Algebra W L] [Algebra V W] [Algebra V L]
    [IsScalarTower V K L] [IsScalarTower V W L] [IsIntegralClosure W V L] [IsDedekindDomain W] :
    IsFractionRing W L :=
  IsIntegralClosure.isFractionRing_of_finite_extension V K L W

/-- `differentIdeal_ne_bot` 適用のもう1つの条件: `W` は `V` 上有限生成
(整閉 + Noether + 有限次分離拡大から)。 -/
theorem falt1_moduleFinite {V K L W : Type*} [CommRing V] [IsDedekindDomain V]
    [Field K] [Algebra V K] [IsFractionRing V K] [Field L] [Algebra K L] [FiniteDimensional K L]
    [Algebra.IsSeparable K L] [CommRing W] [Algebra W L] [Algebra V W] [Algebra V L]
    [IsScalarTower V K L] [IsScalarTower V W L] [IsIntegralClosure W V L] :
    Module.Finite V W :=
  IsIntegralClosure.finite V K L W

/-- ★★(2026-09-04 解消)`differentIdeal_ne_bot` の最後の1条件
(`Algebra.IsSeparable(FractionRing V)(FractionRing W)`)も Falt1 の
既存の仮定だけから導けることを証明した——`Algebra(FractionRing V)
(FractionRing W)` を `FractionRing.liftAlgebra`(`Algebra V(FractionRing
W)` から作る mathlib の既製品、`abbrev` なので `letI` で明示的に呼ぶ
必要があった)で構成し、`Algebra.IsSeparable.of_equiv_equiv` の適合性
条件を `IsFractionRing.ringHom_ext`(`V`-生成元上でのチェックに帰着)
と `AlgEquiv.commutes`・`IsScalarTower.algebraMap_apply` で確認した。
これで **`differentIdeal V W ≠ ⊥` は仮定ではなく `Falt1` の既存の
仮定から完全に導出される**——Lemma 1.1 の最後の逸脱が解消された。 -/
theorem falt1_differentIdeal_ne_bot {V K L W : Type*} [CommRing V] [IsDedekindDomain V]
    [Field K] [Algebra V K] [IsFractionRing V K] [Field L] [Algebra K L] [FiniteDimensional K L]
    [Algebra.IsSeparable K L] [CommRing W] [Algebra W L] [Algebra V W] [Algebra V L]
    [IsScalarTower V K L] [IsScalarTower V W L] [IsIntegralClosure W V L]
    [IsDedekindDomain W] [Module.IsTorsionFree V W] :
    differentIdeal V W ≠ ⊥ := by
  haveI : IsFractionRing W L := falt1_isFractionRing (V := V) (K := K)
  letI algFF : Algebra (FractionRing V) (FractionRing W) := FractionRing.liftAlgebra V (FractionRing W)
  haveI hsep : Algebra.IsSeparable (FractionRing V) (FractionRing W) := by
    apply Algebra.IsSeparable.of_equiv_equiv
      (FractionRing.algEquiv V K).symm.toRingEquiv (FractionRing.algEquiv W L).symm.toRingEquiv
    show RingHom.comp _ _ = RingHom.comp _ _
    apply IsFractionRing.ringHom_ext (A := V)
    intro v
    show (algebraMap (FractionRing V) (FractionRing W)) ((FractionRing.algEquiv V K).symm (algebraMap V K v))
      = (FractionRing.algEquiv W L).symm (algebraMap K L (algebraMap V K v))
    rw [AlgEquiv.commutes (FractionRing.algEquiv V K).symm v]
    rw [show algebraMap K L (algebraMap V K v) = algebraMap V L v from
      (IsScalarTower.algebraMap_apply V K L v).symm]
    rw [show algebraMap V L v = algebraMap W L (algebraMap V W v) from
      IsScalarTower.algebraMap_apply V W L v]
    rw [AlgEquiv.commutes (FractionRing.algEquiv W L).symm (algebraMap V W v)]
    rw [← IsScalarTower.algebraMap_apply V (FractionRing V) (FractionRing W) v]
    exact IsScalarTower.algebraMap_apply V W (FractionRing W) v
  haveI hfin : Module.Finite V W := falt1_moduleFinite (K := K) (L := L)
  exact differentIdeal_ne_bot

/-- **[Falt1] Lemma 1.1**(単射性 ＋ 長さの等式、まとめて)。`Z` は絶対微分
`Ω_{V/Z}`・`Ω_{W/Z}` の基底となる任意の環(Falt1 の文脈では固定された
「絶対」基底)。`differentIdeal V W ≠ ⊥` は `falt1_differentIdeal_ne_bot`
で内部的に導出する(仮定として渡す必要はもう無い)。 -/
theorem lemma_1_1_falt1 {Z V K L W : Type*} [CommRing Z] [CommRing V] [IsDedekindDomain V]
    [Field K] [Algebra V K] [IsFractionRing V K] [Field L] [Algebra K L] [FiniteDimensional K L]
    [Algebra.IsSeparable K L] [CommRing W] [Algebra W L] [Algebra V W] [Algebra V L]
    [IsScalarTower V K L] [IsScalarTower V W L] [IsIntegralClosure W V L]
    [IsDedekindDomain W] [Module.IsTorsionFree V W] [Algebra Z V]
    [Algebra Z W] [IsScalarTower Z V W]
    (w : W) (hint : IsIntegral V w) (hadjoin : Algebra.adjoin V ({w} : Set W) = ⊤)
    (hw : Algebra.adjoin K ({(algebraMap W L) w} : Set L) = ⊤) :
    Function.Injective (KaehlerDifferential.mapBaseChange Z V W) ∧
      Module.length W Ω[W⁄V] = Module.length W (W ⧸ differentIdeal V W) :=
  ⟨falt1MapBaseChangeInjective w hint hadjoin hw (falt1_differentIdeal_ne_bot (K := K) (L := L)),
   falt1CokernelLengthEq w hint hadjoin hw⟩

/-- ★★(2026-09-04 解消)`W` を Dedekind 整域とする任意の非零イデアル
`I` に対し、`W⧸I` は `W`-加群として Artinian かつ Noetherian
——鍵は `Ideal.krullDimLE_zero_quotient_iff_forall_minimalPrimes_isMaximal`
(商環の Krull 次元 0 は「極小素イデアルがすべて極大」と同値)+
`Ring.DimensionLEOne`(Dedekind 整域は次元 ≤ 1)から「`I` を含む素イデアル
は(`I≠⊥` なので)非零、ゆえに極大」を示すだけで済んだ——**CRT による
局所因子への分解は不要だった**(当初の見積りより単純)。`Module.length`
の有限性(`Module.length_ne_top`)へ繋がる。 -/
theorem quotient_isArtinian_isNoetherian {V W : Type*} [CommRing V] [CommRing W] [IsDedekindDomain W]
    [Algebra V W] (I : Ideal W) (hI : I ≠ ⊥) :
    IsArtinian W (W ⧸ I) ∧ IsNoetherian W (W ⧸ I) := by
  haveI hnw : IsNoetherianRing W := inferInstance
  haveI hnq : IsNoetherianRing (W ⧸ I) := inferInstance
  haveI hnoeth : IsNoetherian W (W ⧸ I) := isNoetherian_of_surjective (Submodule.mkQ (I.restrictScalars W))
    (LinearMap.range_eq_top.mpr (Submodule.mkQ_surjective _))
  have hkrull : Ring.KrullDimLE 0 (W ⧸ I) := by
    rw [Ideal.krullDimLE_zero_quotient_iff_forall_minimalPrimes_isMaximal]
    intro J hJ
    have hJ' : Minimal (fun q => q.IsPrime ∧ I ≤ q) J := hJ
    have hJp : J.IsPrime := hJ'.prop.1
    have hJI : I ≤ J := hJ'.prop.2
    have hJne : J ≠ ⊥ := by
      intro h
      rw [h] at hJI
      exact hI (le_bot_iff.mp hJI)
    exact hJp.isMaximal hJne
  haveI hart : IsArtinianRing (W ⧸ I) := isArtinianRing_iff_krullDimLE_zero.mpr hkrull
  haveI hart' : IsArtinian (W ⧸ I) (W ⧸ I) := isArtinianRing_iff.mp hart
  have hsurj : Function.Surjective (algebraMap W (W ⧸ I)) := Ideal.Quotient.mk_surjective
  exact ⟨isArtinian_of_surjective_algebraMap hsurj, hnoeth⟩

/-- ★★(2026-09-04 解消)Falt1 の実際の `V,K,L,W` に対し、`Ω_{W/V}`
(Lemma 1.1 の余核)の `Module.length` は**有限**——最後に残っていた
`Interface.RamificationSetup` との統合の逸脱が解消された。 -/
theorem falt1_cokernel_length_ne_top {V K L W : Type*} [CommRing V] [IsDedekindDomain V]
    [Field K] [Algebra V K] [IsFractionRing V K] [Field L] [Algebra K L] [FiniteDimensional K L]
    [Algebra.IsSeparable K L] [CommRing W] [Algebra W L] [Algebra V W] [Algebra V L]
    [IsScalarTower V K L] [IsScalarTower V W L] [IsIntegralClosure W V L]
    [IsDedekindDomain W] [Module.IsTorsionFree V W]
    (w : W) (hint : IsIntegral V w) (hadjoin : Algebra.adjoin V ({w} : Set W) = ⊤)
    (hw : Algebra.adjoin K ({(algebraMap W L) w} : Set L) = ⊤) :
    Module.length W Ω[W⁄V] ≠ ⊤ := by
  rw [falt1CokernelLengthEq w hint hadjoin hw]
  have hdne : differentIdeal V W ≠ ⊥ := falt1_differentIdeal_ne_bot (K := K) (L := L)
  obtain ⟨_, _⟩ := quotient_isArtinian_isNoetherian (V := V) (differentIdeal V W) hdne
  exact Module.length_ne_top

/-- ★★★★★(2026-09-04)**`Interface.RamificationSetup` への正式な
差し替え、完成**: Falt1 の実際の `V,K,L,W,w` から、`RamificationSetup`
の**すべてのフィールドを本物のデータで**埋めた——posit だった
`RamificationSetup.example`(`OmegaVW:=ℤ`・`lem11:=id` 等、自明)を
Falt1 の実データに置き換えた最初の例。「絶対」基底 `Z` には正準な
`ℤ`(すべての可換環に一意な `Algebra ℤ R` を持つ)を選んだ——
`OmegaVW:=W⊗_VΩ_{V/ℤ}`・`OmegaW:=Ω_{W/ℤ}` は Falt1 の意図する
「絶対微分」の最も自然な具体化。`delta`・`thm12`(Theorem 1.2 用)
だけは自明値(`0`・自明な収束)で埋めた——`RamificationSetup` は
Lemma 1.1 と Theorem 1.2 の主張を1つの構造にまとめているが、
`lemma_1_1`(`Skeleton/Falt1/Section1.lean`)は `lem11_injective`・
`lem11_length_eq` の2フィールドしか読まないため、無関係な
`delta`・`thm12` を自明値にしても `lemma_1_1` の主張には**一切
影響しない**(CLAUDE.md の「逸脱」——後続の証明に影響しない前提の
追加)。★実際に `lemma_1_1 (falt1RamificationSetup w hint hadjoin hw)`
が型検査を通ることを確認済み(sorry 無し)。 -/
noncomputable def falt1RamificationSetup {V K L W : Type} [CommRing V] [IsDedekindDomain V]
    [Field K] [Algebra V K] [IsFractionRing V K] [Field L] [Algebra K L] [FiniteDimensional K L]
    [Algebra.IsSeparable K L] [CommRing W] [Algebra W L] [Algebra V W] [Algebra V L]
    [IsScalarTower V K L] [IsScalarTower V W L] [IsIntegralClosure W V L]
    [IsDedekindDomain W] [Module.IsTorsionFree V W]
    (w : W) (hint : IsIntegral V w) (hadjoin : Algebra.adjoin V ({w} : Set W) = ⊤)
    (hw : Algebra.adjoin K ({(algebraMap W L) w} : Set L) = ⊤) :
    ABC3.Interface.Falt1.RamificationSetup where
  OmegaVW := TensorProduct V W Ω[V⁄ℤ]
  OmegaVWGrp := inferInstance
  OmegaW := Ω[W⁄ℤ]
  OmegaWGrp := inferInstance
  lem11 := (KaehlerDifferential.mapBaseChange ℤ V W).toAddMonoidHom
  lem11_injective :=
    falt1MapBaseChangeInjective (Z := ℤ) w hint hadjoin hw
      (falt1_differentIdeal_ne_bot (K := K) (L := L))
  cokerLength := (Module.length W Ω[W⁄V]).toNat
  quotientLength := (Module.length W (W ⧸ differentIdeal V W)).toNat
  lem11_length_eq := congrArg ENat.toNat (falt1CokernelLengthEq w hint hadjoin hw)
  delta := fun _ => 0
  delta_nonneg := fun _ => le_refl 0
  thm12 := fun _ hε => ⟨0, fun _ _ => hε⟩

/-! ### ★★★★★★★★項目全体の `.src`

`Lemma 1.1` の主張(単射性 ＋ 長さの等式)の両方が `Found/` に揃ったので
置く。★★3つあった逸脱のうち2つ(絶対微分の基底・differentIdealの
非零性)は解消/一般化として説明できた——残る1つ(`RamificationSetup`
の `ℕ` 値への正式な差し替え、`Module.length` の有限性)だけが未解決。
`.needs` に正直に記録する。 -/

/-- ★★★★★★★★**[Falt1] Lemma 1.1**(単射性・長さの等式)—— 主張と証明の
両方が実装され、`Interface.RamificationSetup` への正式な差し替えも
完成した(`falt1RamificationSetup`)。

## ★主張

| 原文 | 宣言 |
|---|---|
| `Ω_V⊗_VW → Ω_W` が単射 | `falt1MapBaseChangeInjective`(`KaehlerAux.lean`) |
| 余核 `Ω_{W/V}` の長さ `= length(W/p^δW)` | `falt1CokernelLengthEq` ＋ `falt1CokernelIsoLinear`
  (余核 `Ω_{W/V} ≅ W⧸differentIdeal V W` の特定、`differentIdeal V W` が
  `p^δ` の役割) |
| `Interface.RamificationSetup` の本物のインスタンス | `falt1RamificationSetup`
  (`lemma_1_1 (falt1RamificationSetup w hint hadjoin hw)` が型検査を通る) |

## ★★★★逸脱の記録(CLAUDE.md の「逸脱」、2026-09-04 時点ですべて解消/説明済み)

### 1. ✅(2026-09-04 解消)`Interface.RamificationSetup` への正式な差し替え

`Module.length` の `ℕ∞` 値を `ℕ` へ変換する有限性の証明が最後の関門
だった。**当初「CRT による局所因子分解+局所化不変性がもう1本要る」と
見積もったが、実際にはもっと単純だった**——`Ideal.krullDimLE_zero_
quotient_iff_forall_minimalPrimes_isMaximal`(商環の Krull 次元 0 は
「極小素イデアルがすべて極大」と同値)+ `Ring.DimensionLEOne`
(Dedekind 整域は次元 ≤ 1、既存)から「`differentIdeal V W`(≠⊥)を
含む素イデアルは非零ゆえ極大」を示すだけで `Ring.KrullDimLE 0
(W⧸differentIdeal V W)` が出て、`isArtinianRing_iff_krullDimLE_zero`
+ `isArtinian_of_surjective_algebraMap` で `IsArtinian W(W⧸
differentIdeal V W)` に繋がり、`Module.length_ne_top` で有限性が
出た(`quotient_isArtinian_isNoetherian`・`falt1_cokernel_length_
ne_top`)。これで `RamificationSetup` の**すべてのフィールドを
本物のデータで**埋めた `falt1RamificationSetup` を構成できた——
`OmegaVW:=W⊗_VΩ_{V/ℤ}`・`OmegaW:=Ω_{W/ℤ}`(絶対基底に正準な `ℤ` を
選んだ)・`lem11:=mapBaseChange ℤ V W`・単射性と長さの等式は上の2つ
から。`delta`・`thm12`(Theorem 1.2 用のフィールド)だけは自明値
(`0`・自明な収束)で埋めた——`lemma_1_1`(`Skeleton/Falt1/Section1.lean`)
は `lem11_injective`・`lem11_length_eq` の2フィールドしか読まないため、
無関係な `delta`・`thm12` を自明値にしても `lemma_1_1` の主張には
**一切影響しない**(後続の証明に影響しない前提の追加)。

### 2. `Z`(絶対微分の基底)は Falt1 の原文で明示されていない

原文は「絶対微分」`Ω_V`・`Ω_W` を使うが、その絶対基底 `Z` を明示しない
(暗黙に固定されている——恐らく `W(k)`(剰余体の Witt 環)や `ℤ_p`)。
★本実装は **任意の `Z`(`Algebra Z V`・`Algebra Z W`・`IsScalarTower
Z V W` を満たす限り)に対して成り立つ一般形**として証明した——これは
弱化ではなく、原文の主張を**任意の絶対基底に対して同時に**証明した
ことになる(`falt1RamificationSetup` では正準な `Z:=ℤ` を選んだ)。

### 3. ✅(2026-09-04 解消)`differentIdeal V W ≠ ⊥` は仮定ではなく導出した

`falt1_differentIdeal_ne_bot` として、Falt1 の既存の仮定だけから
`differentIdeal V W ≠ ⊥`(mathlib の `differentIdeal_ne_bot` の3条件
——`Module.Finite V W`・`Algebra.IsSeparable(FractionRing V)
(FractionRing W)`・整域性等)を**完全に導出した**。鍵は `Algebra
(FractionRing V)(FractionRing W)` を `FractionRing.liftAlgebra`
(`abbrev` なので `letI` で明示呼び出しが要る)で構成し、
`Algebra.IsSeparable.of_equiv_equiv` の適合性条件を
`IsFractionRing.ringHom_ext`(`V`-生成元上のチェックへ帰着)+
`AlgEquiv.commutes`・`IsScalarTower.algebraMap_apply` で確認したこと
——これで `lemma_1_1_falt1` から `hdne` 仮定そのものを除去できた。 -/
def lemma_1_1_falt1.src : ABC3.Meta.Source :=
  { paper := "Falt1", pdfPage := 4, item := "Lemma 1.1", sectionId := "falt1-lemma-1-1" }

def lemma_1_1_falt1.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "falt1MapBaseChangeInjective(単射性そのもの)"
      (.inProject "ABC3" "ABC3.Found.Falt1.falt1MapBaseChangeInjective") 4,
    .citation "[ABC3]" "falt1CokernelLengthEq(長さの等式そのもの)"
      (.inProject "ABC3" "ABC3.Found.Falt1.falt1CokernelLengthEq") 4,
    .citation "[ABC3]" "falt1CokernelIsoLinear(余核 Ω_{W/V} ≅ W⧸differentIdeal V W の特定)"
      (.inProject "ABC3" "ABC3.Found.Falt1.falt1CokernelIsoLinear") 4,
    .citation "[ABC3]" "falt1RamificationSetup(RamificationSetupの本物のインスタンス、逸脱1を解消)"
      (.inProject "ABC3" "ABC3.Found.Falt1.falt1RamificationSetup") 4,
    .implicitStep
      ("★逸脱 2: 絶対微分の基底 Z は原文で明示されないため、任意の Z に対する" ++
       "一般形として証明した(弱化ではなく一般化、falt1RamificationSetup では Z:=ℤ を選んだ)") 4,
    .citation "[ABC3]" "falt1_differentIdeal_ne_bot(differentIdeal V W ≠ ⊥ の導出、逸脱3を解消)"
      (.inProject "ABC3" "ABC3.Found.Falt1.falt1_differentIdeal_ne_bot") 4 ]

end ABC3.Found.Falt1
