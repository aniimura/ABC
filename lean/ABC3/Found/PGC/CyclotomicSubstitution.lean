import ABC3.Found.PGC.ErrorSplitAssembly
import ABC3.Found.PGC.ZetaStepRatio
import ABC3.Found.PGC.ZetaUnitFactor
import ABC3.Found.PGC.ResidueUnitNorm

/-!
# [pGC] 円分塔への代入 —— 3 点のうち 1 点を**定理で供給**、2 点を 1 段まで落とした

## 持ち場（前波の棚卸し表で「★具体層」と書いた 3 行）

`hwn : ‖w‖ = ‖π‖^p` ／ `hvalK`（全分岐）／ `hfr : finrank K M = p`。

## ★結果

| 項目 | 状態 | 中身 |
|---|---|---|
| `hwn` | ★**定理として供給できた**（§2 `norm_w_eq_pow`） | 材料はすべて自分の定理（`ZetaUnitFactor.norm_pow_sub_one_eq` ＋ `ZetaStepRatio.norm_zeta_step`） |
| `hvalK` | ★**1 段下に落ちた**（§3 `valK_of_step`、3 行） | 残るのは「`E₁` の値群が `‖μ‖^ℤ`」＝ `TotallyRamifiedValueGroup.exists_zpow_norm` を `ℚ_p ⊂ E₁` に当てる段 |
| `hfr` | ★**半分は定理**（§4 `le_finrank_of_valK`） | 残るのは `finrank ≤ p` の側だけ（§4 `finrank_eq_of_le`） |

★`hwn` が落ちたのは、★**`w` が「単数倍の差」でしか `μ` と違わない**からである
（`w = ξ^b − 1`、`p ∤ b`）。★`ZetaUnitFactor.norm_pow_sub_one_eq` は前の波で
「使い道が無いかもしれない」と思いながら書いた補題だったが、★ここで効いた。

## ★開いて分かったこと —— 名指しされた 2 ファイルは**使えなかった**

持ち場は `CyclotomicJumpsVerified.lean` / `CyclotomicNumbersVerified.lean` を
「未読の具体層」として挙げていた。★開いた。★**中身は `decide` / `norm_num` レベルの
数値検算**（`(4 : ZMod 27)^9 = 1` など、`p = 3, n = 3` の模型の跳びの数値）であり、
★**一般の円分次数（`finrank`）や値群の定理は 1 つも無い**。⇒ 私の 3 点には使えない。
★「証拠の場所」が外れることもある、という記録として残す。

## ★木の断定の検算（4 例目）

`CyclotomicJumpsVerified.lean` の docstring は
「★おまけ: `ℚ₂(ζ₁₆)/ℚ₂` の break は **`[1, 3, 7]`**、`e = 8`、
`|G_u| = [8,8,4,4,2,2,2,2,1]`」と書いている。
★`tools/zeta-tower-check.py` で**測り直した**:

```
=== L = Q_2(zeta_2^4)   e = deg E = 8 ===
  |G_u| : [8, 8, 4, 4, 2, 2, 2, 2, 1, …]
  ★breaks : [1, 3, 7]
```

★**一致した**（`ℚ₃(ζ₂₇)` の `[2, 8]`、`ℚ₃(ζ₈₁)` の `[2, 8, 26]` も同時に追認）。

## ★次の 1 点に要る部品（測った）

`finrank ≤ p` は「`ζ^p = ξ ∈ E₁` だから `ζ` は `X^p − ξ` の根」から出る:

- `IntermediateField.adjoin.finrank : IsIntegral K x → finrank K K⟮x⟯ = (minpoly K x).natDegree`
  —— ★**在る**（`FieldTheory/IntermediateField/Adjoin/Basic.lean:468`）。
- 次数の上からの評価は ★**`minpoly.min (pmonic) (hp : aeval x p = 0) : degree (minpoly A x) ≤ degree p`**
  （`FieldTheory/Minpoly/Basic.lean:134`）。
- ★**罠**: `minpoly.natDegree_le` は `natDegree ≤ finrank` で ★**向きが逆**。名前が近いので
  引っかかる（#297「形が嘘」の親戚）。

★`IntermediateField` を触るので #59/#69 の領域に入る。★費用は書かない（開いたが試していない）。

## 逸脱の記録

- §3 `valK_of_step` は `K` と `F` の 2 つだけを使う（3 層をまたがない）。★#59 を避けるため、
  「1 段下の値群」を**仮説として受ける**形にした。
- ★`p = 2` は対象外（本鎖の結論が偽になるため）。
-/

namespace ABC3.Found.PGC

namespace CyclotomicSubstitution

/-! ## §1 `‖ζ − 1‖ < 1`（`p` 進で素元は真に小さい） -/

section Small

/-- `‖ζ−1‖^{φ} = ‖p‖ < 1` から `‖ζ−1‖ < 1`。 -/
theorem norm_sub_one_lt_one {F : Type*} [NormedField F] [IsUltrametricDist F]
    [CharZero F] {p m : ℕ} [Fact p.Prime] {ζ : F}
    (hζ : IsPrimitiveRoot ζ (p ^ (m + 1)))
    (hlt : ‖(p : F)‖ < 1) (hne : (p : F) ≠ 0) :
    ‖ζ - 1‖ < 1 := by
  have hpow := ZetaSubOnePrime.norm_zeta_sub_one_pow hζ hlt hne
  by_contra hcon
  have hge : 1 ≤ ‖ζ - 1‖ := not_lt.mp hcon
  have hone : (1 : ℝ) ≤ ‖ζ - 1‖ ^ (Nat.totient (p ^ (m + 1))) := one_le_pow₀ hge
  rw [hpow] at hone
  exact absurd hone (not_le.mpr hlt)

end Small

/-! ## §2 ★`hwn : ‖w‖ = ‖π‖^p` の供給 -/

section W

/-- ★★`w = ξ^b − 1`（`ξ = ζ_{p^{m+1}}`、`p ∤ b`）のノルムは `‖ζ − 1‖^p`。

★これが `MainPartCoeffs` / `ErrorSplitAssembly` の `hwn : ‖w‖ = ‖π‖^p` である。
★材料はすべて自分の定理: `ZetaUnitFactor.norm_pow_sub_one_eq`（単数倍は位を変えない）と
`ZetaStepRatio.norm_zeta_step`（1 段の比はちょうど `p`）。 -/
theorem norm_w_eq_pow {F : Type*} [NormedField F] [IsUltrametricDist F] [CharZero F]
    {p m b : ℕ} [Fact p.Prime] {ζ ξ : F}
    (hζ : IsPrimitiveRoot ζ (p ^ (m + 2))) (hξ : IsPrimitiveRoot ξ (p ^ (m + 1)))
    (hlt : ‖(p : F)‖ < 1) (hne : (p : F) ≠ 0) (hb : ¬ p ∣ b) :
    ‖ξ ^ b - 1‖ = ‖ζ - 1‖ ^ p := by
  have hp : Nat.Prime p := Fact.out
  have hξ1 : ‖ξ‖ = 1 :=
    ZetaUnitFactor.norm_root_of_unity_one (pow_pos hp.pos _) hξ.pow_eq_one
  have hbn : ‖(b : F)‖ = 1 := ResidueUnitNorm.norm_natCast_eq_one_ultra hp hlt hb
  have hsmall : ‖ξ - 1‖ < 1 := norm_sub_one_lt_one hξ hlt hne
  rw [ZetaUnitFactor.norm_pow_sub_one_eq (le_of_eq hξ1) hbn hsmall,
    ← ZetaStepRatio.norm_zeta_step hζ hξ hlt hne]

end W

/-! ## §3 ★`hvalK` を「1 段下の値群」に落とす -/

section ValK

variable {K F : Type*} [Field K] [NormedField F] [Algebra K F]

/-- ★`hvalK`（`Γ_{E₁} ⊆ ‖π‖^{pℤ}`）は、★**`E₁` の値群が `‖μ‖^ℤ`** であることと
`‖μ‖ = ‖π‖^p` から出る（3 行）。

★後者は §2 と同じ `ZetaStepRatio.norm_zeta_step`。
★前者は `TotallyRamifiedValueGroup.exists_zpow_norm` を**1 段下**（`ℚ_p ⊂ E₁`）に
当てたもので、★そこには `finrank ℚ_p E₁ = φ(p^{m+1})` が要る（★未着手）。 -/
theorem valK_of_step {π μ : F} {p : ℕ} (hμ : ‖μ‖ = ‖π‖ ^ p)
    (hstep : ∀ a : K, a ≠ 0 → ∃ k : ℤ, ‖algebraMap K F a‖ = ‖μ‖ ^ k) :
    ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K F a‖ = ‖π‖ ^ ((p : ℤ) * m) := by
  intro a ha
  obtain ⟨k, hk⟩ := hstep a ha
  exact ⟨k, by rw [hk, hμ, ← zpow_natCast, ← zpow_mul]⟩

end ValK

/-! ## §4 ★`hfr : finrank = p` の**半分**は `hvalK` から出る -/

section Finrank

/-- ★`p ≤ finrank K M` は `hvalK` だけから出る（`π^0,…,π^{p−1}` が 1 次独立だから）。

★木の `TotallyRamified.card_le_finrank_of_ne_mod`（`TotallyRamifiedValueGroup.lean:…`）に
代入するだけ。★残る半分（`finrank ≤ p`）は次の 1 点である。 -/
theorem le_finrank_of_valK {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M]
    [Algebra K M] [FiniteDimensional K M] {π : M} {p : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≠ 1)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m)) :
    p ≤ Module.finrank K M := by
  have hcard : Fintype.card (Fin p) ≤ Module.finrank K M := by
    refine TotallyRamified.card_le_finrank_of_ne_mod hπ0 hπ1 hvalK
      (fun l : Fin p => π ^ (l : ℕ)) (fun l => ((l : ℕ) : ℤ)) ?_ ?_
    · intro l; simp [zpow_natCast]
    · intro i j hij
      exact TotallyRamifiedLayer.not_dvd_sub_of_lt i.isLt j.isLt (fun h => hij (Fin.ext h))
  simpa using hcard

/-- ★`hfr` は「`finrank ≤ p`（次数の上からの評価）」だけに落ちる。 -/
theorem finrank_eq_of_le {K M : Type*} [Field K] [NormedField M] [IsUltrametricDist M]
    [Algebra K M] [FiniteDimensional K M] {π : M} {p : ℕ}
    (hπ0 : 0 < ‖π‖) (hπ1 : ‖π‖ ≠ 1)
    (hvalK : ∀ a : K, a ≠ 0 → ∃ m : ℤ, ‖algebraMap K M a‖ = ‖π‖ ^ ((p : ℤ) * m))
    (hle : Module.finrank K M ≤ p) :
    Module.finrank K M = p :=
  le_antisymm hle (le_finrank_of_valK hπ0 hπ1 hvalK)

end Finrank

/-! ## §5 使っている公理の一覧 -/

#print axioms norm_sub_one_lt_one
#print axioms norm_w_eq_pow
#print axioms valK_of_step
#print axioms le_finrank_of_valK
#print axioms finrank_eq_of_le

end CyclotomicSubstitution

end ABC3.Found.PGC
