import ABC3.Found.PGC.LubinTateRamificationBookkeeping
import ABC3.Found.PGC.FixedRingAdjoinIso
import ABC3.Found.PGC.AbelianClosure
import ABC3.Found.PGC.AbelianDecomposition
import ABC3.Found.PGC.AbelianSplitUnramified

/-!
# 局所 Kronecker-Weber の定理 —— Yoshida 2008 Theorem 6.15 の組み立て

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) Theorem 6.15(物理 p.17)。

原文 (Yoshida08 p.17):
> Theorem 6.15. (Local Kronecker-Weber theorem) Every finite abelian extension of a local field K is a Lubin-Tate extension, i.e. K^LT = K^ab.

## 位置づけ —— 道 B の最後の節点

道 B の 6 節点はすべて着地している。本ファイルは**それらを繋ぐだけ**である。

| 節点 | 内容 | ファイル |
|---|---|---|
| B1 | `K^ab` の定義と `K^LT ≤ K^ab` | `AbelianClosure.lean` |
| B2+B3 | `σ` の構成と `K^ab = K^ur·E_σ` | `AbelianClosureSplit.lean` |
| B4 | Prop 6.14(`n = 1`) | `LubinTateRamificationBookkeeping.lean` |
| B5 | `K′ ⊂ K^m_x` | `AbelianSubfieldInLubinTate.lean` |
| 穴 1 | 上付き分岐群の商への降下 | `LubinTateQuotientDescent.lean` |
| 穴 2 | 固定環と整数環の同一視 | `FixedRingAdjoinIso.lean` |
| **B6** | **本ファイル** | |

★★**到達点は `exists_lubinTate_eq_abelianClosure`（§6、仮説なし）**である:

    ∀ K, ∃ E, Normal K E ∧ E ⊓ K^ur = ⊥ ∧ Gal(E/K) ≅ 𝒪_K^× ∧ E ⊔ K^ur = K^ab

★B1 の `exists_lubinTate_le_abelianClosure` は `≤` しか言えなかった。**本ファイルで等号になる。**

## 原典の証明の 5 段（`.txt` 1180–1200 行を直読して回収した）

    Proof. Take a σ ∈ W(K^LT/K) with v(σ) = n > 0, and let L = K_n. Extend σ arbitrarily
    to σ ∈ W(K^ab/K), and let E_σ ⊂ K^ab be its fixed field. Then E_σ ∩ K^ur = L and E_σ/L is
    totally ramified Galois. … Therefore Gal(K^ab/E_σ) ≅ Gal(K^ur E_σ/E_σ), i.e. K^ab = K^ur E_σ.
    Now set x := Art^{-1}_K(σ). Then K^ram_x ⊂ E_σ by Proposition 5.4. As K^LT = K^ur K^ram_x,
    it suffices to show E_σ ⊂ K^ram_x. Let K′/L be any finite Galois extension contained in E_σ.
    It is totally ramified, and Gal(K′/L)_m = {id} for a large m. Then we have
    Gal(K′K^m_x/L)^m = {id} by Proposition 6.14 and Corollary 6.13(ii), hence
    [K′K^m_x : L] | (q−1)q^{m−1} = [K^m_x : L] by Corollary 6.13(iii), thus K′ ⊂ K^m_x.

本ファイルが担当するのは、この 5 段のうち

1. **`K^ab = K^ur·E_σ`**（B2+B3 をそのまま呼ぶ、§2）
2. **`K^ram_x ⊂ E_σ`**（原典の "by Proposition 5.4"。★`n = 1` では
   `σ|_{K_π} = id`（B2）から**直ちに出る** ——§2 の
   `lubinTateClosure_le_fixedField_zpowers`）
3. **`E_σ/K` は Galois**（原典が「Let K′/L be any finite Galois extension contained in E_σ」と
   書くための前提。`E_σ ≤ K^ab` と「交換子を含む部分群は正規」から、§1 の抽象核 3）
4. **無限次の `E_σ` を有限次 Galois 部分拡大で汲み尽くす**（§1 の抽象核 4）
5. **束の計算 `K^ab = K^ur ⊔ E_σ ≤ K^ur ⊔ K_π`**（§1 の抽象核 1・2）

に加えて、5 段目

6. **`K′ ⊂ K^m_x`**（§3 `exists_le_adjoin_psiGenSeq`）—— 穴 2
   `le_of_two_quotients_upperRamification_of_adjoin_pair` に、B4（Prop 6.14）と
   §1b の「for a large m」と原始元定理を渡して閉じる。

である。★★**仮定は 1 つも残っていない。**

★★**§1b は本ファイルで新しく立てた 1 本である**: 原典の "and Gal(K′/L)_m = {id}
for a large m" は**上付き**フィルターで要るのに、木には**下付き**の
`exists_lowerRamificationGroup_eq_bot` しか無かった。Herbrand 関数の辞書
`G^{φ_G(n)} = G_n` と `G^m` の反単調性で乗り換える（`exists_upperRamificationGroup_eq_bot`）。

## 逸脱の記録

1. **`n = 1` に固定した**（決定 D29）。原文は "Take a σ ∈ W(K^LT/K) with v(σ) = n > 0" と
   書いており、`σ` は**我々が選ぶ**ものである。`n = 1`（したがって `L = K_1 = K`）を
   選んでも一般性を失わない。★**選択の固定であって逸脱ではない**が、
   後で `n > 1` の形が要るときのために記録しておく。
2. **`Art^{-1}_K(σ)` を経由しない。** 原文は `x := Art^{-1}_K(σ)` と置いて
   Proposition 5.4 で `K^ram_x ⊂ E_σ` を出すが、`n = 1` では
   B2 が「`K_π` を各点固定し `K^ur` 上 Frobenius になる `σ`」を直接与えるので、
   `K_π ≤ E_σ` は Galois 接続 1 行で出る（§2）。★**原典より短い道**である。
3. **原典の `K′` と `K^m_x` を単項生成に固定した。** 穴 2 の
   `le_of_two_quotients_upperRamification_of_adjoin_pair` が `K(x₀)` / `K(x₁)` の形を
   要求するためである。`K^m_x` は定義から単項生成、`K′` は有限次分離拡大なので
   原始元定理（§1 抽象核 5）で単項生成になる。★**制限になっていない。**
4. **`m`（＝ `k+1`）の選び方は原典と同じではない。** 原文は
   「`Gal(K′/L)_m = {id}` for a large m」でまず `m` を選び、その `m` で `K^m_x` を取る。
   本ファイルも同じ順序だが、`k` は §1b の `exists_upperRamificationGroupAdjoin_eq_bot`
   が返す値をそのまま使う（B4 は**任意の `k`** で成り立つので取り直しが要らない）。

## 退化の自己検査

* ★**`ht`（完全分岐）を落とすと偽**。Corollary 6.13(ii) は
  `K′K′′/K` が完全分岐であることを仮定する。`E_σ` を捨てて
  「2 つの完全分岐アーベル拡大の合成」で代用すると壊れる ——
  `K = ℚ_p` で `ℚ_p(√p)` と `ℚ_p(√(up))`(`u` は非平方単数)はどちらも完全分岐だが、
  その合成は `ℚ_p(√u)`(不分岐)を含む。★**完全分岐性は合成で保たれない。**
  本ファイルはこれを §3 の `isTotallyRamifiedAdjoin_of_le_fixedField_zpowers`
  （`E_σ` の中に居ることから完全分岐が出る）として明示的に供給する。
  ★**`K′` も `K^m_x` も合成体も、すべて `E_σ` の中に居る**ことを使っており、
  `E_σ` を捨てた「2 つの完全分岐拡大」では代用できない。
* ★**`K^ab` は退化していない**: B1 が 4 通りで示している
  （`abelianClosure_ne_bot` / `abelianClosure_selfField_ne_top`(`ℚ_p` のみ) /
  `le_abelianClosure_iff` の両向き / 位相的閉包を取る理由）。§5 でそれを引く。
* ★**結論は等号**である。`≤` は B1 の
  `lubinTateClosure_sup_unramifiedClosure_le_abelianClosure`、
  `≥` が本ファイルの `abelianClosure_le_sup_of_forall_finiteGalois_le`。
  ★片側だけで「示した」と書かない。
* ★**上付き・下付きを取り違えていない**: 本ファイルに出る分岐群は
  §1b の `upperRamificationGroup`（上付き）だけであり、下付き `G_n` は
  在庫 `exists_lowerRamificationGroup_eq_bot` を引くところにしか現れない。
  ★消費側（穴 2）が要求するのも上付きである。
* ★**両辺は `⊥` でも `⊤` でもない**（§5 の 4 本）。とくに `K_π ⊓ K^ur = ⊥` なので、
  `K^ab = K^ur ⊔ K_π` は「片方が他方を含むだけ」の退化した分解ではない。

## 残っている穴

★**無い。** 道 B は本ファイルで閉じる。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-! ## §1 抽象核 —— 分岐・付値・Lubin-Tate の語彙が 1 つも出てこない 5 本

★いずれも mathlib だけで書ける。 -/

section AbstractCore

/-- ★★**抽象核 1**（束）—— `Ω = U ⊔ (Ω ⊓ E)` と `Ω ⊓ E ≤ P` から `Ω ≤ U ⊔ P`。

原典の「`K^ab = K^ur E_σ`」と「`E_σ ⊂ K^ram_x`」を繋ぐ 1 行がこれである。
★**体でも群でもなく、ただの束**でよい。 -/
theorem le_sup_of_sup_inf_eq {α : Type*} [Lattice α] {Ω U E P : α}
    (hsplit : U ⊔ (Ω ⊓ E) = Ω) (hEP : Ω ⊓ E ≤ P) : Ω ≤ U ⊔ P :=
  calc Ω = U ⊔ (Ω ⊓ E) := hsplit.symm
    _ ≤ U ⊔ P := sup_le_sup_left hEP U

/-- ★★**抽象核 2**（束）—— 抽象核 1 に逆向き `U ⊔ P ≤ Ω` を足すと**等号**になる。

★原典の `K^LT = K^ab` は等号である。`hU`（`K^ur ≤ K^ab`）と `hP`（`K_π ≤ K^ab`）は
B1 が与える。 -/
theorem eq_sup_of_sup_inf_eq {α : Type*} [Lattice α] {Ω U E P : α}
    (hsplit : U ⊔ (Ω ⊓ E) = Ω) (hEP : Ω ⊓ E ≤ P) (hU : U ≤ Ω) (hP : P ≤ Ω) :
    Ω = U ⊔ P :=
  le_antisymm (le_sup_of_sup_inf_eq hsplit hEP) (sup_le hU hP)

variable {k L : Type*} [Field k] [Field L] [Algebra k L]

/-- ★**抽象核 3**（群論 + 無限次 Galois 対応）—— **交換子群を固定部分群に含む中間体は
Galois**。

`Gal(L/k)` の部分群 `H` が `⁅Γ,Γ⁆ ≤ H` を満たせば `H ◁ Γ`
（mathlib の `Subgroup.Normal.of_commutator_le`）、それを無限次 Galois 対応
（mathlib の `InfiniteGalois.normal_iff_isGalois`）に渡すだけ。

★原典が「Let K′/L be any finite Galois extension contained in E_σ」と書けるのは、
`E_σ ≤ K^ab` だから `E_σ/K` 自身が（無限次）Galois だからである。 -/
theorem isGalois_of_commutator_le_fixingSubgroup [IsGalois k L] (E : IntermediateField k L)
    (h : commutator (L ≃ₐ[k] L) ≤ E.fixingSubgroup) : IsGalois k E :=
  (InfiniteGalois.normal_iff_isGalois E).mp (Subgroup.Normal.of_commutator_le _ h)

/-- ★★★**抽象核 4**（無限次を有限次で汲み尽くす）——
`E/k` が正規のとき、**`E` に含まれる有限次 Galois 中間体がすべて `F` に入れば `E ≤ F`**。

`x ∈ E` に対し `N := k(x)` の正規閉包（mathlib の
`FiniteGaloisIntermediateField.adjoin k {x}`、`adjoin_val` で
`normalClosure k k(x) L` に等しい）を取ると、`E` が正規なので
`normalClosure_le_iff_of_normal` が `N ≤ E` を与える。

★原典の「Let K′/L be any finite Galois extension contained in E_σ」から
「`E_σ ⊂ K^ram_x`」へ渡る段がちょうどこれである。 -/
theorem le_of_forall_finiteGalois_le [IsGalois k L] (E F : IntermediateField k L) [Normal k ↥E]
    (h : ∀ N : IntermediateField k L, FiniteDimensional k ↥N → IsGalois k ↥N → N ≤ E → N ≤ F) :
    E ≤ F := by
  intro x hx
  have hxE : IntermediateField.adjoin k ({x} : Set L) ≤ E :=
    IntermediateField.adjoin_simple_le_iff.mpr hx
  have hNE : (FiniteGaloisIntermediateField.adjoin k ({x} : Set L)).toIntermediateField ≤ E := by
    rw [FiniteGaloisIntermediateField.adjoin_val]
    exact IntermediateField.normalClosure_le_iff_of_normal.mpr hxE
  have hxN : x ∈ (FiniteGaloisIntermediateField.adjoin k ({x} : Set L)).toIntermediateField :=
    FiniteGaloisIntermediateField.subset_adjoin k ({x} : Set L) rfl
  exact h _ inferInstance inferInstance hNE hxN

/-- ★**抽象核 5**（原始元）—— **有限次分離的な中間体は単項生成**:
`M ≤ L` が `k` 上有限次かつ分離的なら `M = k(y)` なる `y ∈ L` が在る。

mathlib の `Field.exists_primitive_element` は `∃ α : ↥M, k⟮α⟯ = ⊤` という
**`M` の中の**主張なので、`IntermediateField.lift_adjoin_simple` / `lift_top` で
`L` の中に降ろしている（`lean-idioms.md` #59 の定型 (c)：2 層を `lift` に閉じ込める）。

★★本ファイルの結論には使わないが、§3 の合成体 `K′K^m_x` を
B5 が要求する `K(y)` の形にするのに要る。 -/
theorem exists_adjoin_singleton_eq (M : IntermediateField k L)
    [FiniteDimensional k ↥M] [Algebra.IsSeparable k ↥M] :
    ∃ y : L, IntermediateField.adjoin k ({y} : Set L) = M := by
  obtain ⟨α, hα⟩ := Field.exists_primitive_element k ↥M
  exact ⟨(α : L), by
    rw [← IntermediateField.lift_adjoin_simple k M α, hα, IntermediateField.lift_top k M]⟩

/-- ★**抽象核 6**（Galois 接続）—— `σ` が `E` を各点固定するなら `E ⊆ L^σ`。

原典の "Then K^ram_x ⊂ E_σ by Proposition 5.4" にあたる段が、`n = 1` では
**これだけ**になる（B2 が `σ|_{K_π} = id` を与えるので）。 -/
theorem le_fixedField_zpowers_of_mem_fixingSubgroup {E : IntermediateField k L}
    {σ : L ≃ₐ[k] L} (h : σ ∈ E.fixingSubgroup) :
    E ≤ IntermediateField.fixedField (Subgroup.zpowers σ) :=
  (IntermediateField.le_iff_le _ _).mpr (Subgroup.zpowers_le.mpr h)

end AbstractCore

/-! ## §1b 抽象核（一般の DVR と有限群）—— 上付きフィルターは十分大きい `m` で自明

★★原典の "and Gal(K′/L)_m = {id} for a large m" にあたる 1 本。
★Lubin-Tate も p 進体も出てこない（`B` は DVR、`G` は有限群であればよい）。 -/

section EventuallyBot

open IsLocalRing

/-- ★★★**抽象核 7** —— **`G^{k+1} = {1}` なる `k` が在る**（原典の "for a large m"）。

在庫 `exists_lowerRamificationGroup_eq_bot`（下付きフィルターは十分大きい `N` で自明）と
Herbrand 関数の辞書 `G^{φ_G(n)} = G_n`（`upperRamificationGroup_herbrandPhiGroup`）を
繋ぐだけ:

1. `G_N = {1}` なる `N` を取る。
2. `G^{φ_G(N)} = G_N = {1}`（`ramificationGroupReal_eq_of_mem_Ioc` で実数添字を `ℕ` に戻す）。
3. `φ_G(N) ≤ k + 1` なる `k` を取り、`G^m` の反単調性（`upperRamificationGroup_antitone`）。

★★**下付きから上付きへの「for a large m」の乗り換えは木に無かったので、ここで立てた。**
★`k + 1` の形にしてあるのは `ℕ` の切り詰め引き算を書かないため（`lean-idioms.md` #102）
であり、消費側（B5 系）が `((k+1 : ℕ) : ℝ)` を要求する形にそのまま合う。 -/
theorem exists_upperRamificationGroup_eq_bot {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [Fintype G] [SMulCommClass G A B] [FaithfulSMul G B]
    {α : B} (huni : IsLocalRing.maximalIdeal B = Ideal.span {α})
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) :
    ∃ k : ℕ, upperRamificationGroup G α ((k + 1 : ℕ) : ℝ) = ⊥ := by
  obtain ⟨N, hN⟩ := exists_lowerRamificationGroup_eq_bot (A := A) (G := G) hadj
  have hphi : upperRamificationGroup G α (herbrandPhiGroup G α (N : ℝ)) = ⊥ := by
    rw [upperRamificationGroup_herbrandPhiGroup,
      ramificationGroupReal_eq_of_mem_Ioc (A := A) huni hadj N (by linarith) le_rfl, hN N le_rfl]
  obtain ⟨k, hk⟩ := exists_nat_ge (herbrandPhiGroup G α (N : ℝ))
  refine ⟨k, le_antisymm ?_ bot_le⟩
  have hle : herbrandPhiGroup G α (N : ℝ) ≤ ((k + 1 : ℕ) : ℝ) := by push_cast; linarith
  have hmono := upperRamificationGroup_antitone (G := G) huni hle
  rwa [hphi] at hmono

end EventuallyBot

variable {p : ℕ} [Fact p.Prime]

attribute [local instance] isDiscreteValuationRing_adjoinIntegers

/-- ★★**抽象核 7 の具体化** —— 完全分岐な単項拡大 `K(x₁)/K` について
`Gal(K(x₁)/K)^{k+1} = {id}` なる `k` が在る。

`A := 𝒪_K`、`B := 𝒪_{K(x₁)}`。`hadj`（`𝒪_K[α₁] = 𝒪_{K(x₁)}`）は完全分岐 `ht` から
在庫 `adjoin_uniformizer_eq_top_adjoinIntegers` で出る。

★★これが B5 / 穴 2 の `hbot₁`（原典の `Gal(K′/L)^m = {id}`）を埋める。 -/
theorem exists_upperRamificationGroupAdjoin_eq_bot (K : PAdicLocalField p) (x₁ : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x₁} : Set K.closure))]
    [Fintype ((IntermediateField.adjoin K.carrier ({x₁} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x₁} : Set K.closure)))]
    (ht : IsTotallyRamifiedAdjoin K x₁) {α₁ : adjoinIntegers K x₁}
    (huni₁ : IsLocalRing.maximalIdeal (adjoinIntegers K x₁) = Ideal.span {α₁}) :
    ∃ k : ℕ, upperRamificationGroup ((IntermediateField.adjoin K.carrier ({x₁} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x₁} : Set K.closure)))
      α₁ ((k + 1 : ℕ) : ℝ) = ⊥ :=
  exists_upperRamificationGroup_eq_bot (A := 𝒪[K.carrier]) huni₁
    (adjoin_uniformizer_eq_top_adjoinIntegers K x₁ ht huni₁)

section LubinTate

variable (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))

/-! ## §2 具体層 —— `σ` の構成、`E_σ` の Galois 性、`K_π ⊆ E_σ` -/

/-- ★★**原文 "Take a σ ∈ W(K^LT/K) with v(σ) = n > 0"（`n = 1`）**。

`K_π` を各点固定し、`K^ur` 上は算術 Frobenius になる `σ ∈ Γ_K` が在る。
B2 の `exists_arithFrobenius_lift_fixing` に `F := K_π` を代入するだけ
（`K_π ⊓ K^ur = ⊥` は在庫 `lubinTateClosure_inf_unramifiedClosure`）。 -/
theorem exists_arithFrobenius_lift_fixing_lubinTateClosure :
    ∃ σ : K.absGal,
      σ ∈ (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf).fixingSubgroup ∧
      AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
        (↥(unramifiedClosure K)) σ = arithFrobenius K :=
  exists_arithFrobenius_lift_fixing K (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf)
    (normal_lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf)
    (lubinTateClosure_inf_unramifiedClosure K hq hπmax hπne0 f hf0 hf1 hf)

/-- ★★★**原文 "Then K^ram_x ⊂ E_σ by Proposition 5.4"（`n = 1` 版）**。

`σ` が `K_π` を各点固定するので、Galois 接続で `K_π ⊆ K̄^σ`。
★**Artin 写像も Proposition 5.4 も経由しない**（逸脱の記録 2）。 -/
theorem lubinTateClosure_le_fixedField_zpowers {σ : K.absGal}
    (hσfix : σ ∈ (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf).fixingSubgroup) :
    lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf
      ≤ IntermediateField.fixedField (Subgroup.zpowers σ) :=
  le_fixedField_zpowers_of_mem_fixingSubgroup hσfix

/-- ★★**`K_π ⊆ E_σ`**（`E_σ = K^ab ⊓ K̄^σ`）。上の 2 本と B1 の
`lubinTateClosure_le_abelianClosure` を合わせただけ。 -/
theorem lubinTateClosure_le_abelianClosure_inf_fixedField_zpowers {σ : K.absGal}
    (hσfix : σ ∈ (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf).fixingSubgroup) :
    lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf
      ≤ abelianClosure K ⊓ IntermediateField.fixedField (Subgroup.zpowers σ) :=
  le_inf (lubinTateClosure_le_abelianClosure K hq hπmax hπne0 f hf0 hf1 hf)
    (lubinTateClosure_le_fixedField_zpowers K hq hπmax hπne0 f hf0 hf1 hf hσfix)

end LubinTate

/-- ★★**`K^ab` に含まれる中間体は `K` 上（無限次）Galois** —— 抽象核 3 の具体化。

`(K^ab).fixingSubgroup = ⁅Γ_K,Γ_K⁆‾ ⊇ ⁅Γ_K,Γ_K⁆` なので、`A ≤ K^ab` なら
`A.fixingSubgroup` も交換子群を含む。

★原典が「Let K′/L be any finite Galois extension contained in E_σ」と書けるのは
これによる。 -/
theorem isGalois_of_le_abelianClosure (K : PAdicLocalField p)
    (A : IntermediateField K.carrier K.closure) (hA : A ≤ abelianClosure K) :
    IsGalois K.carrier A := by
  haveI := isGalois_closure K
  refine isGalois_of_commutator_le_fixingSubgroup A ?_
  have h1 : (abelianClosure K).fixingSubgroup = (commutator K.absGal).topologicalClosure := by
    rw [abelianClosure_def]
    exact InfiniteGalois.fixingSubgroup_fixedField
      ⟨(commutator K.absGal).topologicalClosure, Subgroup.isClosed_topologicalClosure _⟩
  calc commutator K.absGal ≤ (commutator K.absGal).topologicalClosure :=
        Subgroup.le_topologicalClosure _
    _ = (abelianClosure K).fixingSubgroup := h1.symm
    _ ≤ A.fixingSubgroup := IntermediateField.fixingSubgroup_le hA

/-- ★**`E_σ/K` は Galois** —— 上の系（`E_σ ≤ K^ab` は `inf_le_left`）。 -/
theorem isGalois_abelianClosure_inf_fixedField_zpowers (K : PAdicLocalField p) (σ : K.absGal) :
    IsGalois K.carrier
      ↥(abelianClosure K ⊓ IntermediateField.fixedField (Subgroup.zpowers σ)) :=
  isGalois_of_le_abelianClosure K _ inf_le_left

/-! ## §3 `E_σ` の中の有限次拡大 —— 原典の "It is totally ramified" と `K′ ⊂ K^m_x`

原典の 5 段目
「Let K′/L be any finite Galois extension contained in E_σ. It is totally ramified,
and Gal(K′/L)_m = {id} for a large m. … thus K′ ⊂ K^m_x」
に対応する。★到達点は `exists_le_adjoin_psiGenSeq`（追加仮定ゼロ）。 -/

/-- ★★**原文 "Then E_σ ∩ K^ur = L"（`n = 1` なので `L = K`）の系** ——
`E_σ` に含まれる中間体は `K^ur` と自明にしか交わらない。

B5 の `fixedField_zpowers_inf_unramifiedClosure_eq_bot` に単調性を足しただけ。 -/
theorem inf_unramifiedClosure_eq_bot_of_le_fixedField_zpowers (K : PAdicLocalField p)
    {σ : K.absGal}
    (hσur : AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
      (↥(unramifiedClosure K)) σ = arithFrobenius K)
    {A : IntermediateField K.carrier K.closure}
    (hA : A ≤ IntermediateField.fixedField (Subgroup.zpowers σ)) :
    A ⊓ unramifiedClosure K = ⊥ :=
  le_antisymm
    (calc A ⊓ unramifiedClosure K
        ≤ IntermediateField.fixedField (Subgroup.zpowers σ) ⊓ unramifiedClosure K :=
          inf_le_inf_right _ hA
      _ = ⊥ := fixedField_zpowers_inf_unramifiedClosure_eq_bot K hσur)
    bot_le

/-- ★★★**原文 "It is totally ramified"** —— `E_σ` に含まれる単項拡大は完全分岐。

★★**これを落とすと Corollary 6.13(ii) が偽になる**（モジュール docstring の
`ℚ_p(√p)` / `ℚ_p(√(up))` の反例）。★`E_σ` を捨てられない理由がここにある。 -/
theorem isTotallyRamifiedAdjoin_of_le_fixedField_zpowers (K : PAdicLocalField p)
    {σ : K.absGal}
    (hσur : AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
      (↥(unramifiedClosure K)) σ = arithFrobenius K)
    (y : K.closure)
    (hy : IntermediateField.adjoin K.carrier ({y} : Set K.closure)
      ≤ IntermediateField.fixedField (Subgroup.zpowers σ)) :
    IsTotallyRamifiedAdjoin K y :=
  isTotallyRamifiedAdjoin_of_inf_unramifiedClosure_eq_bot K y
    (inf_unramifiedClosure_eq_bot_of_le_fixedField_zpowers K hσur hy)

/-- ★★**原文 "Gal(K′/L) is abelian"（Corollary 6.13(iii) の入力）** ——
`K^ab` に含まれる正規中間体の Galois 群は可換。B1 の `le_abelianClosure_iff` の
`⇒` 向きそのもの。 -/
theorem forall_mul_comm_of_le_abelianClosure (K : PAdicLocalField p)
    (A : IntermediateField K.carrier K.closure) [Normal K.carrier A]
    (hA : A ≤ abelianClosure K) :
    ∀ a b : (A ≃ₐ[K.carrier] A), a * b = b * a :=
  (le_abelianClosure_iff K A).mp hA

section LubinTate2

variable (K : PAdicLocalField p)
    [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
    [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (hπne0 : π ≠ 0)
    (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
    (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))

/-- ★**`K^m_x ≤ K_π`** —— 原典の "As K^LT = K^ur K^ram_x" の、我々の設定での形。
`K^m_x = K(x_{m+1})` は Lubin-Tate 塔の段そのものである。 -/
theorem adjoin_psiGenSeq_le_lubinTateClosure (m : ℕ) :
    IntermediateField.adjoin K.carrier
        ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure)
      ≤ lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf := by
  rw [lubinTateClosure_eq_adjoin_genSet]
  exact adjoin_psiGenSeq_le_adjoin_genSet K hq hπmax hπne0 f hf0 hf1 hf m

/-- ★★★★**原文 "Let K′/L be any finite Galois extension contained in E_σ … Then we have
Gal(K′K^m_x/L)^m = {id}" の入力を丸ごと供給する**。

`E_σ` に含まれる有限次中間体 `N`（原典の `K′`）と Lubin-Tate 塔の段
`K^m_x = K(x_{m+1})` に対し、その**合成体 `K′K^m_x` は単項生成 `K(y)`** であり、

* **完全分岐**（原文の "It is totally ramified"。★合成体が `E_σ` に入るから出る。
  モジュール docstring の反例のとおり、`E_σ` を捨てるとここが壊れる）
* **Galois 群は可換**（Corollary 6.13(iii) = Hasse-Arf の入力）

★★**これが B5 `le_of_two_quotients_upperRamification_of_finrank` の
`ht` / `habel` / `hsup` をそのまま埋める。** 残るのは `hbot` / `hbot′` / `ϖ` / `ϖ′` /
`hdeg` で、`hdeg` は穴 1 の `finrank_adjoin_iteratedLubinTatePsi_succ` が既に持っている。
★**本ファイルの結論には使わない**（穴 2 が入る次の節点のための材料である）。 -/
theorem exists_adjoin_eq_sup_adjoin_psiGenSeq {σ : K.absGal}
    (hσfix : σ ∈ (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf).fixingSubgroup)
    (hσur : AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
      (↥(unramifiedClosure K)) σ = arithFrobenius K)
    (N : IntermediateField K.carrier K.closure) [FiniteDimensional K.carrier ↥N]
    (hN : N ≤ abelianClosure K ⊓ IntermediateField.fixedField (Subgroup.zpowers σ))
    (m : ℕ) :
    ∃ y : K.closure,
      IntermediateField.adjoin K.carrier ({y} : Set K.closure)
          = N ⊔ IntermediateField.adjoin K.carrier
              ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure) ∧
      IsTotallyRamifiedAdjoin K y ∧
      ∀ a b : (↥(IntermediateField.adjoin K.carrier ({y} : Set K.closure)) ≃ₐ[K.carrier]
        ↥(IntermediateField.adjoin K.carrier ({y} : Set K.closure))), a * b = b * a := by
  haveI := isGalois_closure K
  haveI : FiniteDimensional K.carrier
      ↥(IntermediateField.adjoin K.carrier
        ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure)) :=
    (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).hfd
  have hKmπ : IntermediateField.adjoin K.carrier
      ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure)
      ≤ lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf :=
    adjoin_psiGenSeq_le_lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf m
  have hMσ : (N ⊔ IntermediateField.adjoin K.carrier
      ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure) :
        IntermediateField K.carrier K.closure)
      ≤ IntermediateField.fixedField (Subgroup.zpowers σ) :=
    sup_le (hN.trans inf_le_right)
      (hKmπ.trans (lubinTateClosure_le_fixedField_zpowers K hq hπmax hπne0 f hf0 hf1 hf hσfix))
  have hMab : (N ⊔ IntermediateField.adjoin K.carrier
      ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure) :
        IntermediateField K.carrier K.closure)
      ≤ abelianClosure K :=
    sup_le (hN.trans inf_le_left)
      (hKmπ.trans (lubinTateClosure_le_abelianClosure K hq hπmax hπne0 f hf0 hf1 hf))
  obtain ⟨y, hy⟩ := exists_adjoin_singleton_eq (N ⊔ IntermediateField.adjoin K.carrier
    ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure))
  refine ⟨y, hy, isTotallyRamifiedAdjoin_of_le_fixedField_zpowers K hσur y (hy.trans_le hMσ), ?_⟩
  haveI : IsGalois K.carrier ↥(IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :=
    isGalois_of_le_abelianClosure K _ (hy.trans_le hMab)
  haveI : Normal K.carrier ↥(IntermediateField.adjoin K.carrier ({y} : Set K.closure)) :=
    IsGalois.to_normal
  exact forall_mul_comm_of_le_abelianClosure K _ (hy.trans_le hMab)

def exists_le_adjoin_psiGenSeq.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Theorem 6.15", sectionId := "thm-6-15" }

/-- ★★★★★★★★★★★★**原文 "Let K′/L be any finite Galois extension contained in E_σ …
thus K′ ⊂ K^m_x"** —— 仮定は `K(x₁) ⊆ E_σ` **だけ**。

原文 (Yoshida08 p.17):
> Theorem 6.15. (Local Kronecker-Weber theorem) Every finite abelian extension of a local field K is a Lubin-Tate extension, i.e. K^LT = K^ab.

★★**これが `hfin` を消す本体である。** 段取りは原典の 4 文をそのまま:

| 原文 | 供給元 |
|---|---|
| "It is totally ramified" | §3 `isTotallyRamifiedAdjoin_of_le_fixedField_zpowers`（`E_σ ⊓ K^ur = ⊥`） |
| "Gal(K′/L)_m = {id} for a large m" | §1b `exists_upperRamificationGroupAdjoin_eq_bot` |
| "by Proposition 6.14" | B4 `upperRamificationGroup_torsionGen_eq_bot`（追加仮定ゼロ） |
| "and Corollary 6.13(ii)" | 穴 2 `le_of_two_quotients_upperRamification_of_adjoin_pair` |
| "by Corollary 6.13(iii)" | 同上（`hdeg` は穴 1 `finrank_adjoin_iteratedLubinTatePsi_succ`） |
| "thus K′ ⊂ K^m_x" | 結論 |

合成体 `K′K^m_x` を単項生成 `K(x)` にするのは §1 の抽象核 5（原始元定理）。
`m` は原典の "for a large m" で決まる `k` を使う（`K^m_x = K(x_{k+1})`）——
★**同じ `k` が `hbot₀`（B4、任意の `k` で成立）と `hdeg` にも使える**ので、
`k` の取り直しは要らない。

★★**追加仮定はゼロである**（`htriv` / `hbot` / `hdeg` はすべて外れている）。 -/
theorem exists_le_adjoin_psiGenSeq {σ : K.absGal}
    (hσfix : σ ∈ (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf).fixingSubgroup)
    (hσur : AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
      (↥(unramifiedClosure K)) σ = arithFrobenius K)
    (x₁ : K.closure)
    (hx₁ : IntermediateField.adjoin K.carrier ({x₁} : Set K.closure)
      ≤ abelianClosure K ⊓ IntermediateField.fixedField (Subgroup.zpowers σ)) :
    ∃ m : ℕ, IntermediateField.adjoin K.carrier ({x₁} : Set K.closure)
      ≤ IntermediateField.adjoin K.carrier
        ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf m).pt} : Set K.closure) := by
  haveI := isGalois_closure K
  -- (1) 原典の `K′ = K(x₁)`: 完全分岐で、Galois で、上付きフィルターは `k+1` で自明
  haveI : Fintype ((IntermediateField.adjoin K.carrier ({x₁} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x₁} : Set K.closure))) :=
    Fintype.ofFinite _
  haveI : IsGalois K.carrier ↥(IntermediateField.adjoin K.carrier ({x₁} : Set K.closure)) :=
    isGalois_of_le_abelianClosure K _ (hx₁.trans inf_le_left)
  haveI : Normal K.carrier ↥(IntermediateField.adjoin K.carrier ({x₁} : Set K.closure)) :=
    IsGalois.to_normal
  have ht₁ : IsTotallyRamifiedAdjoin K x₁ :=
    isTotallyRamifiedAdjoin_of_le_fixedField_zpowers K hσur x₁ (hx₁.trans inf_le_right)
  obtain ⟨a₁, ha₁⟩ := IsDiscreteValuationRing.exists_irreducible (adjoinIntegers K x₁)
  have huni₁ : IsLocalRing.maximalIdeal (adjoinIntegers K x₁) = Ideal.span {a₁} :=
    (IsDiscreteValuationRing.irreducible_iff_uniformizer a₁).mp ha₁
  obtain ⟨k, hbot₁⟩ := exists_upperRamificationGroupAdjoin_eq_bot K x₁ ht₁ huni₁
  refine ⟨k, ?_⟩
  -- (2) 原典の `K^m_x = K(x_{k+1})`: B4 が `Gal(K^m_x/K)^{k+1} = {id}` を与える
  have hψ₀ := (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).hψ
  have htor₀ := (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).hn
  have hmem₀ := (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).hmem
  haveI := (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).hfd
  haveI : Fintype ((IntermediateField.adjoin K.carrier
        ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).pt} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier
        ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).pt} : Set K.closure))) :=
    Fintype.ofFinite _
  have hx₀π := adjoin_psiGenSeq_le_lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf k
  have hx₀ab := hx₀π.trans (lubinTateClosure_le_abelianClosure K hq hπmax hπne0 f hf0 hf1 hf)
  have hx₀σ := hx₀π.trans
    (lubinTateClosure_le_fixedField_zpowers K hq hπmax hπne0 f hf0 hf1 hf hσfix)
  haveI : IsGalois K.carrier ↥(IntermediateField.adjoin K.carrier
      ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).pt} : Set K.closure)) :=
    isGalois_of_le_abelianClosure K _ hx₀ab
  haveI : Normal K.carrier ↥(IntermediateField.adjoin K.carrier
      ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).pt} : Set K.closure)) := IsGalois.to_normal
  have hbot₀ := upperRamificationGroup_torsionGen_eq_bot K hq hπmax hπne0 f hf0 hf1 hf k
    (Nat.le_add_left 1 k) _ hψ₀ htor₀ hmem₀
  have huni₀ := (IsDiscreteValuationRing.irreducible_iff_uniformizer
      (torsionGen K hq hπmax hπne0 f hf0 hf1 hf (k + 1) _ htor₀ hmem₀)).mp
    (irreducible_torsionGen K hq hπmax hπne0 f hf0 hf1 hf (k + 1) (Nat.le_add_left 1 k) _
      hψ₀ htor₀ hmem₀)
  have hdeg := finrank_adjoin_iteratedLubinTatePsi_succ K hq hπmax hπne0 f hf0 hf1 hf k _
    hψ₀ htor₀ hmem₀
  -- (3) 合成体 `K′K^m_x` を単項生成 `K(x)` にする（抽象核 5＝原始元定理）
  obtain ⟨x, hx⟩ := exists_adjoin_singleton_eq
    (IntermediateField.adjoin K.carrier
        ({(psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).pt} : Set K.closure)
      ⊔ IntermediateField.adjoin K.carrier ({x₁} : Set K.closure))
  have hsup := hx.symm
  have h₁ := le_sup_left.trans_eq hsup
  have h₂ := le_sup_right.trans_eq hsup
  have hxσ := hx.trans_le (sup_le hx₀σ (hx₁.trans inf_le_right))
  have hxab := hx.trans_le (sup_le hx₀ab (hx₁.trans inf_le_left))
  have ht : IsTotallyRamifiedAdjoin K x :=
    isTotallyRamifiedAdjoin_of_le_fixedField_zpowers K hσur x hxσ
  haveI : IsGalois K.carrier ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    isGalois_of_le_abelianClosure K _ hxab
  haveI : Normal K.carrier ↥(IntermediateField.adjoin K.carrier ({x} : Set K.closure)) :=
    IsGalois.to_normal
  haveI : Fintype ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) :=
    Fintype.ofFinite _
  obtain ⟨a, ha⟩ := IsDiscreteValuationRing.exists_irreducible (adjoinIntegers K x)
  have huni := (IsDiscreteValuationRing.irreducible_iff_uniformizer a).mp ha
  have habel := forall_mul_comm_of_le_abelianClosure K _ hxab
  haveI : Fintype (((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      ⧸ (IntermediateField.restrict h₁).fixingSubgroup) := Fintype.ofFinite _
  haveI : Fintype (((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
      ⧸ (IntermediateField.restrict h₂).fixingSubgroup) := Fintype.ofFinite _
  letI := quotientMulSemiringActionOfTrivial
    (↥(fixedRing ↥(adjoinIntegers K x) (IntermediateField.restrict h₁).fixingSubgroup))
    (smul_fixedRing_eq_self (B := ↥(adjoinIntegers K x))
      (H := (IntermediateField.restrict h₁).fixingSubgroup))
  letI := quotientMulSemiringActionOfTrivial
    (↥(fixedRing ↥(adjoinIntegers K x) (IntermediateField.restrict h₂).fixingSubgroup))
    (smul_fixedRing_eq_self (B := ↥(adjoinIntegers K x))
      (H := (IntermediateField.restrict h₂).fixingSubgroup))
  exact le_of_two_quotients_upperRamification_of_adjoin_pair K x
    (psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf k).pt ht huni habel x₁ h₁ h₂ hsup
    (fun ρ c => quotientSMul_mk_fixedRing ρ c) (fun ρ c => quotientSMul_mk_fixedRing ρ c)
    k huni₀ hbot₀ huni₁ hbot₁ hdeg

/-! ## §4 主定理 -/

def abelianClosure_le_sup_of_forall_finiteGalois_le.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Theorem 6.15", sectionId := "thm-6-15" }

/-- ★★★★★★★★**Yoshida 2008 Theorem 6.15 の `≥` 側**: `K^ab ≤ K_π ⊔ K^ur = K^LT`。

原文 (Yoshida08 p.17):
> Theorem 6.15. (Local Kronecker-Weber theorem) Every finite abelian extension of a local field K is a Lubin-Tate extension, i.e. K^LT = K^ab.

★★**`≤` 側は B1（`lubinTateClosure_sup_unramifiedClosure_le_abelianClosure`）が
既に持っている。本補題が出すのは `≥` 側である。取り違えないこと。**

段取り（原典の 5 段のうち 4 段）:
1. B2+B3 の `fixedField_commutator_eq_unramifiedClosure_sup_inf_fixedField_zpowers` で
   `K^ab = K^ur ⊔ E_σ`。
2. `E_σ ≤ K^ab` と抽象核 3 で `E_σ/K` は Galois（したがって正規）。
3. 抽象核 4 で `E_σ ≤ K_π` を「`E_σ` に含まれる有限次 Galois 中間体がすべて `K_π` に入る」
   （仮定 `hfin`）に帰着。
4. 抽象核 1（束）で `K^ab ≤ K^ur ⊔ K_π`。

★逸脱 3: `hfin` は仮定である（原典の "thus K′ ⊂ K^m_x" = B5）。
穴 2 が入れば B5 の `le_of_two_quotients_upperRamification_of_finrank` が埋める。 -/
theorem abelianClosure_le_sup_of_forall_finiteGalois_le {σ : K.absGal}
    (hσur : AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
      (↥(unramifiedClosure K)) σ = arithFrobenius K)
    (hfin : ∀ N : IntermediateField K.carrier K.closure,
      FiniteDimensional K.carrier ↥N → IsGalois K.carrier ↥N →
      N ≤ abelianClosure K ⊓ IntermediateField.fixedField (Subgroup.zpowers σ) →
      N ≤ lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf) :
    abelianClosure K
      ≤ (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
          IntermediateField K.carrier K.closure) := by
  haveI := isGalois_closure K
  haveI := isGalois_abelianClosure_inf_fixedField_zpowers K σ
  haveI : Normal K.carrier
      ↥(abelianClosure K ⊓ IntermediateField.fixedField (Subgroup.zpowers σ)) :=
    IsGalois.to_normal
  have hsplit : unramifiedClosure K ⊔
      (abelianClosure K ⊓ IntermediateField.fixedField (Subgroup.zpowers σ))
      = abelianClosure K :=
    fixedField_commutator_eq_unramifiedClosure_sup_inf_fixedField_zpowers K hσur
  have hEP : abelianClosure K ⊓ IntermediateField.fixedField (Subgroup.zpowers σ)
      ≤ lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf :=
    le_of_forall_finiteGalois_le _ _ hfin
  have := le_sup_of_sup_inf_eq hsplit hEP
  rwa [sup_comm] at this

def abelianClosure_eq_sup_of_forall_finiteGalois_le.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Theorem 6.15", sectionId := "thm-6-15" }

/-- ★★★★★★★★★★★★**Theorem 6.15 の等号形（`E_σ ⊆ K_π` を仮定に置いた段階）**。

原文 (Yoshida08 p.17):
> Theorem 6.15. (Local Kronecker-Weber theorem) Every finite abelian extension of a local field K is a Lubin-Tate extension, i.e. K^LT = K^ab.

`≤` は B1、`≥` は上の `abelianClosure_le_sup_of_forall_finiteGalois_le`。
★**両側そろって初めて等号になる。**
★下の `abelianClosure_eq_lubinTateClosure_sup_unramifiedClosure` で `hfin` は外れる。 -/
theorem abelianClosure_eq_sup_of_forall_finiteGalois_le {σ : K.absGal}
    (hσur : AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
      (↥(unramifiedClosure K)) σ = arithFrobenius K)
    (hfin : ∀ N : IntermediateField K.carrier K.closure,
      FiniteDimensional K.carrier ↥N → IsGalois K.carrier ↥N →
      N ≤ abelianClosure K ⊓ IntermediateField.fixedField (Subgroup.zpowers σ) →
      N ≤ lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf) :
    abelianClosure K
      = (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
          IntermediateField K.carrier K.closure) :=
  le_antisymm
    (abelianClosure_le_sup_of_forall_finiteGalois_le K hq hπmax hπne0 f hf0 hf1 hf hσur hfin)
    (lubinTateClosure_sup_unramifiedClosure_le_abelianClosure K hq hπmax hπne0 f hf0 hf1 hf)

def forall_finiteGalois_le_lubinTateClosure.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Theorem 6.15", sectionId := "thm-6-15" }

/-- ★★★★★★★★★★**原文 "it suffices to show E_σ ⊂ K^ram_x" の中身** ——
`E_σ` に含まれる有限次 Galois 中間体はすべて `K_π` に入る。

原文 (Yoshida08 p.17):
> Theorem 6.15. (Local Kronecker-Weber theorem) Every finite abelian extension of a local field K is a Lubin-Tate extension, i.e. K^LT = K^ab.

`N` を原始元定理（§1 抽象核 5）で `K(x₁)` に書き直し、§3 の
`exists_le_adjoin_psiGenSeq` に渡すだけ。★★**追加仮定はゼロ。** -/
theorem forall_finiteGalois_le_lubinTateClosure {σ : K.absGal}
    (hσfix : σ ∈ (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf).fixingSubgroup)
    (hσur : AlgEquiv.restrictNormalHom (F := K.carrier) (K₁ := K.closure)
      (↥(unramifiedClosure K)) σ = arithFrobenius K) :
    ∀ N : IntermediateField K.carrier K.closure,
      FiniteDimensional K.carrier ↥N → IsGalois K.carrier ↥N →
      N ≤ abelianClosure K ⊓ IntermediateField.fixedField (Subgroup.zpowers σ) →
      N ≤ lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf := by
  haveI := isGalois_closure K
  intro N hfd _ hN
  haveI := hfd
  obtain ⟨x₁, hx₁⟩ := exists_adjoin_singleton_eq N
  obtain ⟨m, hm⟩ := exists_le_adjoin_psiGenSeq K hq hπmax hπne0 f hf0 hf1 hf hσfix hσur x₁
    (hx₁.trans_le hN)
  exact hx₁.symm.trans_le
    (hm.trans (adjoin_psiGenSeq_le_lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf m))

def abelianClosure_eq_lubinTateClosure_sup_unramifiedClosure.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Theorem 6.15", sectionId := "thm-6-15" }

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**Yoshida 2008 Theorem 6.15（局所 Kronecker-Weber の定理）—— `K^LT = K^ab`。
★★追加仮定ゼロ。**

原文 (Yoshida08 p.17):
> Theorem 6.15. (Local Kronecker-Weber theorem) Every finite abelian extension of a local field K is a Lubin-Tate extension, i.e. K^LT = K^ab.

★★★**これが道 B の到達点である。**

* `≤`（`K^LT ≤ K^ab`）は B1 `lubinTateClosure_sup_unramifiedClosure_le_abelianClosure`。
* `≥`（`K^ab ≤ K^LT`）が本ファイル。`σ` は B2 が構成し（"Take a σ"）、
  `K^ab = K^ur·E_σ` は B3、`E_σ ⊆ K_π` は §3 の `exists_le_adjoin_psiGenSeq`。

★仮定に残っているのは **Lubin-Tate の設定そのもの**（`hq`・`hπmax`・`hπne0`・
`f` が Lubin-Tate 級数であること）だけであり、それも下の
`exists_lubinTate_eq_abelianClosure` で消える。 -/
theorem abelianClosure_eq_lubinTateClosure_sup_unramifiedClosure :
    abelianClosure K
      = (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
          IntermediateField K.carrier K.closure) := by
  obtain ⟨σ, hσfix, hσur⟩ :=
    exists_arithFrobenius_lift_fixing_lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf
  exact abelianClosure_eq_sup_of_forall_finiteGalois_le K hq hπmax hπne0 f hf0 hf1 hf hσur
    (forall_finiteGalois_le_lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf hσfix hσur)

/-! ## §5 退化の自己検査 -/

/-- ★**自己検査 1** —— 結論の右辺 `K^LT = K_π · K^ur` は `⊥` に潰れていない。

★これがないと「`K^ab = K^LT`」は**両辺 `⊥` の自明な等式**でありうる。
`K^ur ≠ K`（B1 の `unramifiedClosure_ne_bot`、`Gal(K^ur/K) ≅ Ẑ` が非自明だから）から。
★左辺 `K^ab ≠ ⊥` は B1 の `abelianClosure_ne_bot` が既に持っている。 -/
theorem lubinTateClosure_sup_unramifiedClosure_ne_bot :
    (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) ≠ ⊥ := fun h =>
  unramifiedClosure_ne_bot K (le_antisymm (le_trans le_sup_right (le_of_eq h)) bot_le)

/-- ★**自己検査 2** —— 結論の両辺は `⊤`（＝ `K̄`）にも潰れていない。

`K^ab ≠ K̄` は `Γ_K` が非可換であることに同値で、B1 は `K = ℚ_p` について
`abelianClosure_selfField_ne_top` で示している。★一般の `K` についての
`Γ_K` 非可換はまだ木に無いので、`hne` を仮定として受け取る形にしてある。 -/
theorem lubinTateClosure_sup_unramifiedClosure_ne_top
    (hne : abelianClosure K ≠ ⊤)
    (heq : abelianClosure K
      = (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
          IntermediateField K.carrier K.closure)) :
    (lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊔ unramifiedClosure K :
      IntermediateField K.carrier K.closure) ≠ ⊤ := heq ▸ hne

/-- ★**自己検査 4** —— `K_π` と `K^ur` は**自明にしか交わらない**（在庫の再掲）。

★これで「`K^ab = K^ur ⊔ K_π`」が「片方が他方を含むだけ」の退化した分解では
ないことが分かる（`K_π ≠ ⊥` かつ `K^ur ≠ ⊥` かつ `K_π ⊓ K^ur = ⊥`）。 -/
theorem lubinTateClosure_inf_unramifiedClosure_eq_bot :
    lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf ⊓ unramifiedClosure K = ⊥ :=
  lubinTateClosure_inf_unramifiedClosure K hq hπmax hπne0 f hf0 hf1 hf

end LubinTate2

/-! ## §6 仮説なしの形 —— 原典の Theorem 6.15 そのもの -/

def exists_lubinTate_eq_abelianClosure.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 17, item := "Theorem 6.15", sectionId := "thm-6-15" }

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**Yoshida 2008 Theorem 6.15（局所 Kronecker-Weber の定理）—— 仮説なしの形。**

原文 (Yoshida08 p.17):
> Theorem 6.15. (Local Kronecker-Weber theorem) Every finite abelian extension of a local field K is a Lubin-Tate extension, i.e. K^LT = K^ab.

任意の p 進局所体 `K` について、`K̄` の中に `K` 上正規な `K_π` が在って

* `K_π ⊓ K^ur = K`
* `Gal(K_π/K) ≅ 𝒪_K^×`
* ★★**`K_π · K^ur = K^ab`（等号）**

★★★**B1 の `exists_lubinTate_le_abelianClosure` は `≤` しか言えなかった。
本定理はそれを等号にする。** 道 B（B1–B6 + 穴 1 + 穴 2）はこれで閉じる。

★`Λ_∞` の選び方（`π`・`f`）に依らない形にするため、`K_π` は `∃` の内側に
閉じ込めてある（`exists_lubinTate_le_abelianClosure` と同じ流儀）。
★Lubin-Tate の設定（adic 完備性・剰余体の有限性と位数・素元・Lubin-Tate 級数）は
すべて本リポジトリで既に構築済みなので、`AbelianDecomposition.lean` の
`exists_lubinTateUnramified_decomposition` と同じ組み上げで消える。 -/
theorem exists_lubinTate_eq_abelianClosure (K : PAdicLocalField p) :
    ∃ E : IntermediateField K.carrier K.closure,
      Normal K.carrier E ∧
      E ⊓ unramifiedClosure K = ⊥ ∧
      Nonempty ((E ≃ₐ[K.carrier] E) ≃* (𝒪[K.carrier])ˣ) ∧
      (E ⊔ unramifiedClosure K : IntermediateField K.carrier K.closure) = abelianClosure K := by
  haveI := isAdicComplete_valuationRing K
  haveI := valuationRing_isDVR K
  obtain ⟨π, hπirr⟩ := IsDiscreteValuationRing.exists_irreducible (𝒪[K.carrier])
  have hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π} :=
    (IsDiscreteValuationRing.irreducible_iff_uniformizer π).mp hπirr
  have hπne0 : π ≠ 0 := hπirr.ne_zero
  have hq : Fintype.card 𝓀[K.carrier] = p ^ (absoluteInertiaDegree K) := by
    rw [← Nat.card_eq_fintype_card]
    exact residueCard_eq_pow K
  obtain ⟨f, hf0, hf1, hf⟩ := exists_lubinTateSeries (A := 𝒪[K.carrier]) hq hπmax
  exact ⟨lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf,
    normal_lubinTateClosure K hq hπmax hπne0 f hf0 hf1 hf,
    lubinTateClosure_inf_unramifiedClosure K hq hπmax hπne0 f hf0 hf1 hf,
    ⟨lubinTateClosureGalEquivUnits K hq hπmax hπne0 f hf0 hf1 hf⟩,
    (abelianClosure_eq_lubinTateClosure_sup_unramifiedClosure K hq hπmax hπne0 f hf0 hf1
      hf).symm⟩

end ABC3.Found.PGC
