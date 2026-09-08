import ABC3.Found.PGC.WildBreakLowerBound
import ABC3.Found.PGC.RamificationJumpBound

/-!
# [pGC] 暴分岐の跳びの**上界** `(p−1)i ≤ p·e` —— これで `k = 0` の出口から分岐の仮説が消えた

`Found/PGC/WildBreakLowerBound.lean` が下界 `1 ≤ t` を閉じ、残った 1 点が上界だった。
★本ファイルはそれを閉じ、★★**`k = 0`(次数 `p` の層)の出口から跳びについての仮説を
1 本残らず消す**(`..._of_uniformizer_deg_p_free`)。

## ★★前波(＝自分)の記述の訂正(名指し)

第 1105 の報告と `WildBreakLowerBound.lean` 冒頭はこう書いた:

> ★残り 1 点は**上界** `(p−1)t ≤ p·e`。その筋は「`f'(π) = ∏_{j≠0}(π − g^jπ)` と
> Eisenstein 多項式の係数評価」で、★下界と違って**最小多項式の係数が整**であることが要る。

★★**「係数が整であることが要る」は誤りである。**
`Found/PGC/RamificationJumpBound.lean:316` の `norm_natCast_le_pow_of_splits` が要求するのは

    Monic / Separable / natDegree = n / aeval π f = 0 / Splits / hbreak(共役が等距離)

だけで、★**係数の整性(Eisenstein 性)も `n` の素数性も一度も使わない**。
同ファイル冒頭 docstring の「★`differentIdeal` を通らず、Eisenstein 性も `n` の素数性も
要らず monic だけで足りた」の方が正しかった。
★元の docstring は他が読んでいるので直さず、ここに訂正として書く。

★**測ったコマンド**: `sed -n '316,345p' lean/ABC3/Found/PGC/RamificationJumpBound.lean`
（`hf : f.Monic` / `hdeg` / `hsep : f.Separable` / `hroot` / `hsp : Splits` / `hbreak` の 6 本のみ）。

## 何を足したか(＝在庫に無かったのはどこか)

`RamificationJumpBound` が要求する `hbreak`

    ∀ a, (minpoly K π).map(algebraMap).IsRoot a → a ≠ π → ‖π − a‖ = ‖π‖^{i+1}

を供給するのに 2 つ要った。★どちらも在庫に無かった。

1. §1 ★**共役はすべて等距離**(`norm_pow_apply_sub_eq`)。
   `‖g^j π − π‖ ≤ ‖g π − π‖` は超距離＋群だけで出る。`p` 素数・`p ∤ j` なら `g^j` も生成元
   なので**両向き**に使えて等号になる。★分岐・付値の語彙が 1 語も出ない。
2. §2 ★**根はすべて `g` の軌道の中**(`exists_pow_of_isRoot`)。
   `IsConjRoot.exists_algEquiv`(`Normal` が要る)＋
   `Nat.card ⟨g⟩ = orderOf g = p = [M:K] = Nat.card Gal(M/K)` ⇒ `⟨g⟩ = ⊤`。
   ★**ノルムが 1 語も出ない**。

## ★#323 の薬が効いた(実測)

§2 の 2 本(`natDegree_minpoly_of_adjoin_eq_top` / `exists_pow_of_isRoot`)は
★**`[Field M]` だけで書いた**(`NormedField` を使わない)。#323 の
「`NormedField M` のまま `Algebra.adjoin` の所属を書くと `isDefEq` で焼き切れる」を踏まないため。
★`IntermediateField` は `K⟮π⟯` **1 つだけ**作り、その場で `⊤` に潰す(#59 の 2 層をまたがない)。

## 出口の仮説の推移(実測)

| 出口 | 跳び `t` についての仮説 |
|---|---|
| `TotallyRamifiedLayer.…_of_uniformizer_deg_p` | `1 ≤ t` かつ `(p−1)t ≤ p·e` |
| `WildBreak.…_of_uniformizer_deg_p_upper`(第 1105) | `(p−1)t ≤ p·e` だけ |
| ★★`WildBreakUpper.…_of_uniformizer_deg_p_free`(本ファイル) | ★**無し** |

★最終形に残るのは
`orderOf g = p` / `[M:K] = p` / `IsGalois K M` / `hiso`(等長) / `hvalK`(全分岐) /
`hnormp`(正規化) / `hπlt`(素元) の **7 本**で、★分岐の**数値**についての条件は 0 本。

## ★閉じていないもの(正確に)

* `IsGalois K M` を仮説に置いた。★`orderOf g = p` と `[M:K] = p` から自動で出るか
  (`Fintype.card Gal ≤ finrank` が一般の有限次拡大で成り立つか)は★**測っていない**。
  ★応用先(`M^P/M^Q`)は Galois なので、この仮説は無料である。
* `hiso`(等長)は `PureStepSetup.lean:283 norm_algEquiv_eq` が供給する(第 1102 で測った)。
* ★**不分岐側**(`f = p`)と★**体の側の翻訳**(`M^P/M^Q` を `IntermediateField` に)は未着手。

## ★在庫の測定(コマンドを残す)

```
#check @IsConjRoot.exists_algEquiv   → [Normal K L] (h : IsConjRoot K x y) : ∃ σ, σ y = x
#check @isConjRoot_iff_aeval_eq_zero → IsIntegral K x → (IsConjRoot K x y ↔ aeval y (minpoly K x) = 0)
#check @Normal.splits / @IsGalois.card_aut_eq_finrank / @Subgroup.eq_top_of_card_eq
#check @mem_powers_iff_mem_zpowers  → [Finite G] : y ∈ powers x ↔ y ∈ zpowers x
#check @IntermediateField.adjoin.finrank → finrank K ↥K⟮x⟯ = (minpoly K x).natDegree
#check @IntermediateField.adjoin_eq_top_of_algebra → Algebra.adjoin F S = ⊤ → adjoin F S = ⊤
★`ZMod.natCast_zmod_eq_zero_iff_dvd` は**無い**。正しくは `ZMod.natCast_eq_zero_iff`。
★`le_or_lt` は Unknown identifier(`by_cases` に書き換えた)。
★`IntermediateField.finrank_top` は `finrank ↥⊤ E = 1` で**向きが違う**。
  `finrank K ↥⊤ = finrank K M` は `IntermediateField.topEquiv` ＋ `LinearEquiv.finrank_eq`。
```
-/

namespace ABC3.Found.PGC

namespace WildBreakUpper

open Finset WildBreak

/-! ## §1 抽象核 —— 共役はすべて等距離 -/

section Equidistant

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★**抽象核** —— `g` が等長なら `‖g^j π − π‖ ≤ ‖g π − π‖`(超距離＋群の性質だけ)。 -/
theorem norm_pow_apply_sub_le (g : M ≃ₐ[K] M) (hiso : ∀ z : M, ‖g z‖ = ‖z‖) (π : M) :
    ∀ j : ℕ, ‖(g ^ j) π - π‖ ≤ ‖g π - π‖ := by
  intro j
  induction j with
  | zero => simp
  | succ j ih =>
      have hrw : (g ^ (j + 1)) π - π = g ((g ^ j) π - π) + (g π - π) := by
        rw [map_sub, ← AlgEquiv.mul_apply, ← pow_succ']
        ring
      rw [hrw]
      refine (IsUltrametricDist.norm_add_le_max _ _).trans (max_le ?_ le_rfl)
      rw [hiso]
      exact ih

/-- ★★**抽象核** —— `orderOf g = p`(素数)で `p ∤ j` なら `‖g^j π − π‖ = ‖g π − π‖`。

★`g^j` も生成元なので上の不等式を**両向き**に使える。★分岐・付値の語彙は 1 語も出ない。 -/
theorem norm_pow_apply_sub_eq {p : ℕ} (hp : p.Prime) (g : M ≃ₐ[K] M) (hg : orderOf g = p)
    (hiso : ∀ z : M, ‖g z‖ = ‖z‖) (π : M) {j : ℕ} (hj : ¬ (p ∣ j)) :
    ‖(g ^ j) π - π‖ = ‖g π - π‖ := by
  haveI : Fact p.Prime := ⟨hp⟩
  refine le_antisymm (norm_pow_apply_sub_le g hiso π j) ?_
  -- `j` の法 `p` 逆元 `k` を取る
  have hjne : ((j : ZMod p)) ≠ 0 := by
    rw [Ne, ZMod.natCast_eq_zero_iff]
    exact hj
  set k : ℕ := ((j : ZMod p)⁻¹).val with hkdef
  have hmod : j * k ≡ 1 [MOD p] := by
    rw [← ZMod.natCast_eq_natCast_iff]
    push_cast
    rw [hkdef, ZMod.natCast_val, ZMod.cast_id]
    field_simp
  have hgjk : (g ^ j) ^ k = g := by
    rw [← pow_mul, ← pow_one g]
    exact (pow_eq_pow_iff_modEq).mpr (by simpa [hg] using hmod)
  have hisoj : ∀ z : M, ‖(g ^ j) z‖ = ‖z‖ := norm_pow_apply g hiso j
  have := norm_pow_apply_sub_le (g ^ j) hisoj π k
  rwa [hgjk] at this

end Equidistant

/-! ## §2 Galois の段 —— ★**ノルムが 1 語も出ない**(#323 の薬: `[Field M]` だけで書く) -/

section GaloisPart

/-- `M = K(π)` なら `deg minpoly = [M:K]`。★`NormedField` を**使わない**形で書く(#323)。 -/
theorem natDegree_minpoly_of_adjoin_eq_top {K M : Type*} [Field K] [Field M] [Algebra K M]
    [FiniteDimensional K M] {π : M} (htop : Algebra.adjoin K ({π} : Set M) = ⊤) :
    (minpoly K π).natDegree = Module.finrank K M := by
  have hint : IsIntegral K π := IsIntegral.of_finite K π
  have htop' : IntermediateField.adjoin K ({π} : Set M) = ⊤ :=
    IntermediateField.adjoin_eq_top_of_algebra K _ htop
  have h1 : Module.finrank K (IntermediateField.adjoin K ({π} : Set M))
      = (minpoly K π).natDegree := IntermediateField.adjoin.finrank hint
  have h2 : Module.finrank K (IntermediateField.adjoin K ({π} : Set M))
      = Module.finrank K M := by
    rw [htop']
    exact LinearEquiv.finrank_eq (IntermediateField.topEquiv (F := K) (E := M)).toLinearEquiv
  rw [← h1, h2]

/-- ★★**`minpoly K π` の根はすべて `g` の軌道の中**。

`IsConjRoot.exists_algEquiv`(`Normal` が要る)で共役を与える `σ` を取り、
`Nat.card ⟨g⟩ = orderOf g = p = [M:K] = Nat.card Gal(M/K)` から `⟨g⟩ = ⊤` を出す。
★ノルムが 1 語も出ない(#323 の薬)。 -/
theorem exists_pow_of_isRoot {K M : Type*} [Field K] [Field M] [Algebra K M]
    [FiniteDimensional K M] [IsGalois K M] {p : ℕ} {π : M} (g : M ≃ₐ[K] M)
    (hg : orderOf g = p) (hnK : Module.finrank K M = p)
    {a : M} (haroot : Polynomial.aeval a (minpoly K π) = 0) :
    ∃ j : ℕ, (g ^ j) π = a := by
  have hint : IsIntegral K π := IsIntegral.of_finite K π
  have hconj : IsConjRoot K π a := (isConjRoot_iff_aeval_eq_zero hint).mpr haroot
  obtain ⟨σ, hσ⟩ := IsConjRoot.exists_algEquiv hconj
  have hcard : Nat.card (Subgroup.zpowers g) = Nat.card (M ≃ₐ[K] M) := by
    rw [Nat.card_zpowers, hg, IsGalois.card_aut_eq_finrank, hnK]
  have htop : Subgroup.zpowers g = ⊤ := Subgroup.eq_top_of_card_eq _ hcard
  have hmem : σ⁻¹ ∈ Subgroup.zpowers g := by rw [htop]; trivial
  obtain ⟨j, hj⟩ := mem_powers_iff_mem_zpowers.mpr hmem
  have hj' : g ^ j = σ⁻¹ := hj
  refine ⟨j, ?_⟩
  rw [hj', ← hσ]
  simp

end GaloisPart
/-! ## §3 ★★★★上界 `(p−1)·i ≤ p·e` の供給 -/

section Upper

variable {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M] [Algebra K M]

/-- ★★★★★**跳びの上界** —— 次数 `p` の全分岐巡回拡大で `(p−1)·i ≤ p·e`。

`RamificationJumpBound.norm_natCast_le_pow_of_splits`(monic だけで足りる。★Eisenstein 性も
係数の整性も要らない)に、§1(共役は等距離)と §2(根は軌道の中)を差し込む。

★★**前波(＝自分)の記述の訂正**: 第 1105 の報告は
「上界には★**最小多項式の係数が整**であることが要る」と書いたが、★**これは誤りである。**
`norm_natCast_le_pow_of_splits` が要求するのは `Monic` / `Separable` / `Splits` /
`natDegree = n` だけで、係数の整性(Eisenstein 性)は**一度も使わない**
(`RamificationJumpBound.lean` 冒頭 docstring の「monic だけで足りた」が正しい)。 -/
theorem sub_one_mul_le_of_totallyRamified [FiniteDimensional K M] [IsGalois K M]
    {p e i : ℕ} (hp : p.Prime) {π : M}
    (g : M ≃ₐ[K] M) (hg : orderOf g = p) (hiso : ∀ z : M, ‖g z‖ = ‖z‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hnK : Module.finrank K M = p)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (heM : ‖(p : M)‖ = ‖π‖ ^ (p * e))
    (hbr : ‖g π - π‖ = ‖π‖ ^ (i + 1)) :
    (p - 1) * i ≤ p * e := by
  have hint : IsIntegral K π := IsIntegral.of_finite K π
  have htop : Algebra.adjoin K ({π} : Set M) = ⊤ :=
    TotallyRamifiedLayer.adjoin_eq_top_of_valK hπ0 (ne_of_lt hπ1) hnK hvalK
  have hdeg : (minpoly K π).natDegree = p := by
    rw [natDegree_minpoly_of_adjoin_eq_top htop, hnK]
  have hgp : g ^ p = 1 := by rw [← hg]; exact pow_orderOf_eq_one g
  have hbreak : ∀ a : M, ((minpoly K π).map (algebraMap K M)).IsRoot a → a ≠ π →
      ‖π - a‖ = ‖π‖ ^ (i + 1) := by
    intro a haroot hane
    have haeval : Polynomial.aeval a (minpoly K π) = 0 := by
      rw [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def] at haroot
      exact haroot
    obtain ⟨j, hj⟩ := exists_pow_of_isRoot (p := p) g hg hnK haeval
    have hnd : ¬ (p ∣ j) := by
      intro ⟨c, hc⟩
      apply hane
      rw [← hj, hc, pow_mul, hgp, one_pow]
      rfl
    rw [norm_sub_rev, ← hj, norm_pow_apply_sub_eq hp g hg hiso π hnd, hbr]
  have hle : ‖(p : M)‖ ≤ ‖π‖ ^ ((p - 1) * i) :=
    norm_natCast_le_pow_of_splits hp.pos hπ0 hπ1 hvalK (minpoly.monic hint) hdeg
      (Algebra.IsSeparable.isSeparable K π) (minpoly.aeval K π)
      (Normal.splits (IsGalois.to_normal) π) hbreak
  exact sub_one_mul_le_of_norm_natCast_eq_pow hπ0 hπ1 heM hle

/-- ★★★★★★**`hupper` の供給** —— 出口が要求する `ℤ` の形。

`t ≤ 0` の場合は `(p−1)t ≤ 0 ≤ p·e` で自明。`1 ≤ t` の場合に上を使う。 -/
theorem hupper_of_totallyRamified [FiniteDimensional K M] [IsGalois K M]
    {p : ℕ} (hp : p.Prime) {π : M}
    (g : M ≃ₐ[K] M) (hg : orderOf g = p) (hiso : ∀ z : M, ‖g z‖ = ‖z‖)
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ < 1) (hnK : Module.finrank K M = p)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (e : ℕ) (t : ℤ) (he : 0 < e) (heM : ‖(p : M)‖ = ‖π‖ ^ (p * e))
    (hbr : ‖g π - π‖ = ‖π‖ ^ (t + 1)) :
    ((p : ℤ) - 1) * t ≤ (p : ℤ) * (e : ℤ) := by
  have hp1 : (1 : ℤ) ≤ (p : ℤ) := by exact_mod_cast hp.one_le
  have hpe : (0 : ℤ) ≤ (p : ℤ) * (e : ℤ) := by positivity
  by_cases ht : 0 < t
  · have hi : ((t.toNat : ℕ) : ℤ) = t := Int.toNat_of_nonneg (le_of_lt ht)
    have hbrN : ‖g π - π‖ = ‖π‖ ^ (t.toNat + 1) := by
      rw [hbr, ← zpow_natCast ‖π‖ (t.toNat + 1)]
      congr 1
      push_cast [hi]
      ring
    have hnat := sub_one_mul_le_of_totallyRamified hp g hg hiso hπ0 hπ1 hnK hvalK heM hbrN
    have hcast : (((p - 1) * t.toNat : ℕ) : ℤ) ≤ ((p * e : ℕ) : ℤ) := by exact_mod_cast hnat
    push_cast [Nat.cast_sub hp.one_le, hi] at hcast
    linarith
  · rw [not_lt] at ht
    nlinarith

end Upper

/-! ## §4 ★★★★★★★出口 —— **分岐についての仮説が 1 本も無い** -/

section Exit

open GainedTowerModel GainedTowerModel.GaloisTower WildBreak TotallyRamifiedLayer

variable {M : Type*} [NormedField M] [IsUltrametricDist M]

/-- ★★★★★★★**`k = 0`(次数 `p` の層)の出口から、分岐についての仮説が全部消えた形**。

残る仮説は
「`g` が位数 `p`」「`[M:K] = p`」「`M/K` が Galois」「`g` が等長」
「全分岐(`hvalK`)」「ノルムの正規化(`hnormp`)」「`π` が素元(`hπlt`)」だけで、
★**跳び `t` についての条件は 1 本も無い**(下界は `WildBreakLowerBound` §4、
上界は本ファイル §3 が定理として供給する)。 -/
theorem exists_norm_sub_algebraMap_le_axDecay_of_uniformizer_deg_p_free
    {K : Type*} [Field K] [Algebra K M] [FiniteDimensional K M] [IsGalois K M]
    {p : ℕ} [Fact p.Prime] {π : M}
    (g : M ≃ₐ[K] M) (hg : orderOf g = p) (hiso : ∀ z : M, ‖g z‖ = ‖z‖)
    (hnK : Module.finrank K M = p)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hnormp : ‖(p : M)‖ = ((p : ℝ))⁻¹) (hπlt : ‖π‖ < 1) (x : M) :
    ∃ y : twr g p 0, ‖x - algebraMap (twr g p 0) M y‖ ≤ axDecay p 1 * ‖g x - x‖ := by
  have hp : p.Prime := Fact.out
  have hπ0 : 0 < ‖π‖ := norm_pos_of_valK (p := p) (n := p) hnormp hvalK
  exact exists_norm_sub_algebraMap_le_axDecay_of_uniformizer_deg_p_upper g hg hiso
    (fun e t he heM hbr =>
      hupper_of_totallyRamified hp g hg hiso hπ0 hπlt hnK hvalK e t he heM hbr)
    hnK hvalK hnormp hπlt x

end Exit


/-! ## §5 `.src`(原典の対応箇所) -/

def norm_pow_apply_sub_eq.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def sub_one_mul_le_of_totallyRamified.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

def exists_norm_sub_algebraMap_le_axDecay_of_uniformizer_deg_p_free.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 6, item := "Corollary 3.1", sectionId := "cor-3-1" }

/-! ## §6 使っている公理の一覧 -/

#print axioms norm_pow_apply_sub_le
#print axioms norm_pow_apply_sub_eq
#print axioms natDegree_minpoly_of_adjoin_eq_top
#print axioms exists_pow_of_isRoot
#print axioms sub_one_mul_le_of_totallyRamified
#print axioms hupper_of_totallyRamified
#print axioms exists_norm_sub_algebraMap_le_axDecay_of_uniformizer_deg_p_free

end WildBreakUpper

end ABC3.Found.PGC
