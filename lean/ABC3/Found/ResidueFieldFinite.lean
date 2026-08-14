import Mathlib.Topology.Algebra.Valued.LocallyCompact
import Mathlib.FieldTheory.Finite.Basic

/-!
# 完備・proper な非アルキメデスノルム体の剰余体は有限、かつ位数は素数冪

論文にも我々のモデルにも固有でない、**一般の**結果。mathlib へ出せる形で書く。

## 経緯(記録として残す)

当初は mathlib の
`Valued.integer.properSpace_iff_..._finite_residueField` を直接使おうとして、
**ノルムのダイヤモンド**で止まった——その補題の `ProperSpace` は
`Valued.toNormedField`(付値から作ったノルム)に関するもので、
こちらの `ProperSpace` はスペクトルノルムに関するもの。
両者は命題としては等しいが定義的に等しくない。

round-trip 補題(`Valued.toNormedField (NormedField.toValued) ≍ 元のノルム`)を
書きにいったが、**それは不要だった**——`CompactSpace 𝒪[K]` を経由すれば、
問題はノルムでなく**位相**の一致に落ちる。そして `NormedField.toValued` は
`toUniformSpace` に元のノルムのものをそのまま使うので、位相は**定義的に**一致する。

教訓: ダイヤモンドは「等しいことを証明する」より「**触れずに済ませる**」方が安い。

## 剰余体の標数(2026-08-14 追加)

有限性だけでは `Interface` の `ResidueCardinality` は埋まらない——位数が **p の冪**
であることが要る。そのために `CharP 𝓀[K] p` を作る。仮説は
**`‖(p : K)‖ < 1`** ひとつ(これが「p は 𝒪[K] の単元でない」を与える)。

### ★先行する測定の訂正(2026-08-14)

`Found/PGC/LocalFieldNorm.lean` の docstring に
「mathlib に `charP_of_prime_eq_zero` に相当する直接の補題は見当たらず(実測)、
`ringChar` 経由で数行書く必要がある」と書いてあったが、**これは誤り**だった。

```
Mathlib/Algebra/CharP/Basic.lean:103
  CharP.charP_iff_prime_eq_zero [Nontrivial R] {p : ℕ} (hp : p.Prime) :
    CharP R p ↔ (p : R) = 0
```

が存在する(同 100 行に `CharP.ringChar_of_prime_eq_zero` もある)。
前回は `Mathlib/Algebra/CharP/Defs.lean` の `ringChar` 名前空間しか見ておらず、
`Basic.lean` を見ていなかった。**「無い」という測定は、探索範囲を書かないと再現できない**
——`LeanStatus.absent` を記録するときの教訓として残す。

したがって `ringChar` を手で回す必要は無く、`charP_iff_prime_eq_zero` 一発で済む。
残る仕事は `(p : 𝓀[K]) = 0`、すなわち **p が `𝒪[K]` の単元でない**ことだけ。
-/

open scoped NormedField Valued

variable {K : Type*} [NontriviallyNormedField K] [IsUltrametricDist K] [ProperSpace K]

omit [ProperSpace K] in
/-- 付値環は閉単位球にほかならない。 -/
theorem integer_eq_closedBall : (𝒪[K] : Set K) = Metric.closedBall 0 1 := by
  ext x
  simp [Valued.integer.mem_iff]

/-- proper なら付値環はコンパクト(閉単位球だから)。 -/
instance compactSpace_integer : CompactSpace 𝒪[K] :=
  isCompact_iff_compactSpace.mp (integer_eq_closedBall (K := K) ▸ isCompact_closedBall 0 1)

/-- **剰余体は有限**。 -/
theorem finite_residueField : Finite 𝓀[K] :=
  letI : (Valued.v : Valuation K NNReal).RankOne :=
    inferInstanceAs (NormedField.valuation (K := K)).RankOne
  (Valued.integer.compactSpace_iff_completeSpace_and_isDiscreteValuationRing_and_finite_residueField
    (K := K) (Γ₀ := NNReal)).mp inferInstance |>.2.2

/-! ## 剰余体の標数と位数 -/

omit [ProperSpace K] in
/-- `𝒪[K]` の元としてのノルムが 1 未満なら、それは `𝒪[K]` の単元でない。

付値の言葉に落として mathlib の
`Valuation.Integer.not_isUnit_iff_valuation_lt_one` を使う。 -/
theorem not_isUnit_of_norm_lt_one {x : 𝒪[K]} (hx : ‖(x : K)‖ < 1) : ¬ IsUnit x := by
  rw [Valuation.Integer.not_isUnit_iff_valuation_lt_one]
  simpa [← NNReal.coe_lt_coe] using hx

omit [ProperSpace K] in
/-- ノルムが 1 未満の自然数は、剰余体では `0`。 -/
theorem natCast_residueField_eq_zero {n : ℕ} (hn : ‖(n : K)‖ < 1) : (n : 𝓀[K]) = 0 := by
  have hcoe : (((n : 𝒪[K]) : K)) = (n : K) := by push_cast; rfl
  rw [← map_natCast (IsLocalRing.residue 𝒪[K]) n, IsLocalRing.residue_eq_zero_iff,
    IsLocalRing.mem_maximalIdeal, mem_nonunits_iff]
  exact not_isUnit_of_norm_lt_one (by rw [hcoe]; exact hn)

omit [ProperSpace K] in
/-- **剰余体の標数**。`‖(n : K)‖ < 1` を満たす素数 `n` は剰余体の標数である。 -/
theorem charP_residueField {n : ℕ} (hn : n.Prime) (h : ‖(n : K)‖ < 1) : CharP 𝓀[K] n :=
  (CharP.charP_iff_prime_eq_zero hn).mpr (natCast_residueField_eq_zero h)

omit [ProperSpace K] in
/-- **剰余体を別の局所環へ移す橋**。付値環が局所環 `A` と環同型なら、剰余体の位数は等しい。

付値環を「ノルムが 1 以下の元の集合」として同定できたとき、そこから先はこれ 1 本で済む
——`A` 側で剰余体が既知なら(例: `A = ℤ_[p]` に対する `PadicInt.residueField`)、
`𝓀[K]` の位数がそのまま出る。

★**ノルムの一致そのものは要らない**。要るのは「付値環が集合として一致すること」だけ。
`Found/PGC/QpResidueField.lean` はこれを使って `ℚ_[p]` の場合を片付けている。 -/
theorem card_residueField_eq_of_ringEquiv {A : Type*} [CommRing A] [IsLocalRing A]
    (e : 𝒪[K] ≃+* A) : Nat.card 𝓀[K] = Nat.card (IsLocalRing.ResidueField A) :=
  Nat.card_congr (IsLocalRing.ResidueField.mapEquiv e).toEquiv

/-- **剰余体の位数は `n` の正の冪**(原文 [pGC] p.3 の「k is the field of q = p^f elements」)。

`0 < f` を落とすと `q = 1`(自明体)を許してしまうので、ここは主張の一部。 -/
theorem card_residueField_eq_prime_pow {n : ℕ} (hn : n.Prime) (h : ‖(n : K)‖ < 1) :
    ∃ f : ℕ, 0 < f ∧ Nat.card 𝓀[K] = n ^ f := by
  haveI : CharP 𝓀[K] n := charP_residueField hn h
  haveI : Finite 𝓀[K] := finite_residueField
  haveI : Fintype 𝓀[K] := Fintype.ofFinite _
  obtain ⟨f, -, hf⟩ := FiniteField.card 𝓀[K] n
  exact ⟨(f : ℕ), f.2, by rw [Nat.card_eq_fintype_card]; exact hf⟩
