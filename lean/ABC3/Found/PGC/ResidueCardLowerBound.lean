import ABC3.Found.PGC.KummerDuality
import ABC3.Found.PGC.UnitsPowPrimeToP
import ABC3.Found.PGC.UnramifiedRootsOfUnity

/-!
# `n ∣ q − 1` のとき `N_n(Γ_K) = n²`(経路 C の下界、G1)

[pGC] Proposition 1.2 への**経路 C**(`ResearchPaper/pgc-goal.md`)は、絶対 Galois 群
`Γ_K` の位相群としての同型類だけから剰余体の位数 `q` を読み取ることを目指す。その
不変量が

  `N_n(Γ) := #Hom_cont(Γ, ℤ/n)`   (`Found/PGC/ContinuousHomCount.lean::contHomCard`)

であり、経路 C の主張 (C-q) は「`p ∤ n` のとき `N_n(Γ_K) = n · gcd(n, q−1)`」、
したがって

  `q − 1 = max { n : p ∤ n, N_n(Γ_K) = n² }`

である。本ファイルはこの `max` の**下界の側**——`n ∣ q−1` なる `n` は
実際に `N_n(Γ_K) = n²` を満たす——を `sorry` 無しで確定させる。
上界(`n ∤ q−1` なら `N_n < n²`)は経路 C のノード F(`InertiaKummer.lean`)側の仕事で、
ここでは扱わない。

## 証明(在庫 3 つの合成。新しい数学は `μ_n ⊆ K` の一行だけ)

1. **`μ_n ⊆ K`**(`exists_isPrimitiveRoot_mem_carrier`)。
   `𝒪_K^×` の `p` と素な捩れ部分群 `μ^{(p')}(𝒪_K)`
   (`Found/PGC/PrimeToPTorsion.lean::primeToPTorsion`)は位数 `q−1` の**巡回**群
   (`card_primeToPTorsion` + `isCyclic_primeToPTorsion`)なので、`n ∣ q−1` から
   位数ちょうど `n` の元が取れる
   (`Found/PGC/UnramifiedRootsOfUnity.lean::exists_orderOf_eq_of_dvd_card`)。
   単射準同型 `μ^{(p')}(𝒪_K) ↪ K^×` で位数は変わらないので、その像は原始 `n` 乗根。
2. **Kummer 双対**(`Found/PGC/KummerDuality.lean::card_contHom_eq_index_powRange`)を
   `F := K.carrier` に当てて
   `N_n(Γ_K) = [K^× : (K^×)^n]`。
3. **B3**(`Found/PGC/UnitsPowPrimeToP.lean::index_powRange_carrierUnits`)で
   `[K^× : (K^×)^n] = n · gcd(n, q−1)`。`n ∣ q−1` だから `gcd(n, q−1) = n`。

★配管について。`Γ_K = K.closure ≃ₐ[K.carrier] K.closure` で `K.closure` は
`AlgebraicClosure K.carrier` の `abbrev`(`Skeleton/PGC/Setup.lean`)なので、
Kummer 双対が述べている `AlgebraicClosure K.carrier ≃ₐ[K.carrier] AlgebraicClosure K.carrier`
との同一視は **`rfl` で済む**(位相・群の instance も一致する)。
`Found/PGC/AdjoinFieldClosure.lean` の `autCongrContinuousMulEquiv` や
`contHomCard_congr` による移送は**不要だった**。

## ★退化の自己検査——落とした仮定は主張を偽にするか

★**この定理で落とせる仮定は `hdvd` ただ 1 つである。** `hn` と `hn0` は `hdvd` から
**従う**ので、「落とす」という操作が意味を持たない。これは持ち場を受け取った時点の
見込み(「`hn` を落とすと B3 が破れて偽になる」)とは違っていたので、そのまま記録する。

* **`hdvd : n ∣ q−1` を落とすと偽**。正しい値は経路 C の (C-q)
  `N_n(Γ_K) = n · gcd(n, q−1)` であって、`n ∤ q−1` なら `gcd(n, q−1) < n` だから
  `N_n < n²` になる。例えば `K = ℚ_5`(`q = 5`, `q − 1 = 4`)・`n = 3` では
  `gcd(3, 4) = 1` だから `N_3(Γ_{ℚ_5}) = 3 ≠ 9`。実際このとき `μ_3 ⊄ ℚ_5` なので、
  上の 1. からして成り立たない。
  ★この不等号こそが経路 C の主張の実体である(`q−1` が `max` として決まる)。
* **`hn : ¬ p ∣ n` は `hdvd` から従う**(`not_dvd_of_dvd_residueCard_sub_one`)。
  `p ∤ q−1`(`PrimeToPTorsion.lean::not_dvd_residueCard_sub_one`)なので、
  `p ∣ n` と `n ∣ q−1` は両立しない。
  ——B3(`index_powRange_carrierUnits`)**単体**では `¬ p ∣ n` は本質的で、落とすと
  偽になる(`n = p` のとき `[K^× : (K^×)^p] = p^{[K:ℚ_p]+1} · #μ_p(K) ≠ p`。実例は
  第 1016 の docstring)。だがそれは B3 の退化であって、本定理の退化ではない。
  本定理では `hdvd` が `hn` を含んでしまうため、`hn` を外した反例は作れない。
* **`hn0 : n ≠ 0` も `hdvd` から従う**(`ne_zero_of_dvd_residueCard_sub_one`)。
  `n = 0` なら `q − 1 = 0` だが `q ≥ 2`(`Fintype.one_lt_card`)なので矛盾。
* 冗長な 2 つを署名に残しているのは、下流(ノード G の総組み立て)が
  `hn` を手元に持ったまま呼べるようにするためである。中身は
  `card_contHom_eq_index_powRange`(`hn0` を要求)と
  `index_powRange_carrierUnits`(`hn` を要求)にそのまま渡している。

## 逸脱の記録

原典 [pGC] は Proposition 1.2 の論拠を Serre の局所類体論(相互律)に投げているが、
経路 C はそれを経由しない(`ResearchPaper/pgc-goal.md` で記録済みの逸脱)。
本ファイルもその方針のまま、Hensel・Teichmüller・素手の Kummer 指標だけで済ませている。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC ABC3.Found
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

/-! ## `μ_n ⊆ K`(`n ∣ q−1` のとき) -/

/-- `p` と素な捩れ部分群 `μ^{(p')}(𝒪_K)` から `K^×` への準同型
(`𝒪_K ↪ K` を単数群に持ち上げただけ)。 -/
noncomputable def primeToPTorsionToUnits (K : PAdicLocalField p) :
    primeToPTorsion K →* (K.carrier)ˣ :=
  (Units.map (algebraMap 𝒪[K.carrier] K.carrier).toMonoidHom).comp (primeToPTorsion K).subtype

theorem primeToPTorsionToUnits_injective (K : PAdicLocalField p) :
    Function.Injective (primeToPTorsionToUnits K) :=
  (Units.map_injective (fun _ _ h => Subtype.ext h)).comp Subtype.val_injective

/-- **★★★★★★★★`n ∣ q − 1` なら `μ_n ⊆ K`**——原始 `n` 乗根が `K` 自身の中にある。

`μ^{(p')}(𝒪_K)` は位数 `q−1` の巡回群なので、その約数 `n` ごとに位数ちょうど `n` の
元がある(`exists_orderOf_eq_of_dvd_card`)。単射準同型 `μ^{(p')}(𝒪_K) ↪ K^×` は
位数を保つので、像は `K` の原始 `n` 乗根である。

`Found/PGC/UnramifiedRootsOfUnity.lean::exists_isPrimitiveRoot_mem_unramifiedClosure`
(`p ∤ m` なら `μ_m ⊆ K^ur`)の**`K` 自身の版**にあたる。あちらは `q^k − 1` を作るために
不分岐拡大を張る必要があったが、こちらは `n ∣ q−1` という仮定のおかげで拡大が要らない。 -/
theorem exists_isPrimitiveRoot_mem_carrier (K : PAdicLocalField p) {n : ℕ}
    (hdvd : n ∣ Nat.card 𝓀[K.carrier] - 1) :
    ∃ ζ : K.carrier, IsPrimitiveRoot ζ n := by
  haveI := isCyclic_primeToPTorsion K
  haveI : Finite (primeToPTorsion K) := Finite.of_equiv _ (primeToPTorsionEquiv K).symm.toEquiv
  obtain ⟨u, hu⟩ := exists_orderOf_eq_of_dvd_card
    (G := primeToPTorsion K) (m := n) (by rw [card_primeToPTorsion]; exact hdvd)
  refine ⟨((primeToPTorsionToUnits K u : (K.carrier)ˣ) : K.carrier), ?_⟩
  rw [IsPrimitiveRoot.coe_units_iff]
  have h : orderOf (primeToPTorsionToUnits K u) = n := by
    rw [orderOf_injective _ (primeToPTorsionToUnits_injective K) u, hu]
  exact h ▸ IsPrimitiveRoot.orderOf _

/-! ## 退化の自己検査——`hn` と `hn0` は `hdvd` から従う

本ファイル冒頭に書いた「落とせる仮定は `hdvd` だけ」を、散文でなく Lean で確定させる。 -/

/-- `n ∣ q − 1` なら `n ≠ 0`(`q ≥ 2` なので `q − 1 ≥ 1`)。 -/
theorem ne_zero_of_dvd_residueCard_sub_one (K : PAdicLocalField p) {n : ℕ}
    (hdvd : n ∣ Nat.card 𝓀[K.carrier] - 1) : n ≠ 0 := by
  rintro rfl
  have h1 : 1 < Fintype.card 𝓀[K.carrier] := Fintype.one_lt_card
  have h2 : Nat.card 𝓀[K.carrier] - 1 = 0 := Nat.eq_zero_of_zero_dvd hdvd
  rw [Nat.card_eq_fintype_card] at h2
  omega

/-- `n ∣ q − 1` なら `p ∤ n`(`p ∤ q − 1` だから)。 -/
theorem not_dvd_of_dvd_residueCard_sub_one (K : PAdicLocalField p) {n : ℕ}
    (hdvd : n ∣ Nat.card 𝓀[K.carrier] - 1) : ¬ p ∣ n := fun hp =>
  not_dvd_residueCard_sub_one K (by rw [← Nat.card_eq_fintype_card]; exact hp.trans hdvd)

/-! ## G1——`N_n(Γ_K) = n²` -/

/-- **★★★★★★★★★★★★★★★★経路 C の下界(G1):`n ∣ q − 1` なら
`N_n(Γ_K) = n²`**。

`Γ_K` からの連続指標の個数が `n²` に**達する**ことを言う。経路 C はこの等式が
`n ∣ q−1` を特徴づける(`q − 1 = max { n : p ∤ n, N_n(Γ_K) = n² }`)ことで
`q` を群論的に回復する。ここはその `max` が実際に取られること——下界——である。

証明は在庫 3 つの合成:`μ_n ⊆ K`(`exists_isPrimitiveRoot_mem_carrier`)⟹
Kummer 双対(`card_contHom_eq_index_powRange`)で `N_n = [K^× : (K^×)^n]` ⟹
B3(`index_powRange_carrierUnits`)で `= n · gcd(n, q−1) = n · n`。

★仮定の非自明性は本ファイル冒頭の「退化の自己検査」に書いた。落とせるのは `hdvd`
だけで、落とすと偽になる(`K = ℚ_5`, `n = 3` で `N_3 = 3 ≠ 9`)。`hn` と `hn0` は
`hdvd` から従う(`not_dvd_of_dvd_residueCard_sub_one` /
`ne_zero_of_dvd_residueCard_sub_one`)ので、外した反例は存在しない。 -/
theorem contHomCard_absGal_of_dvd (K : PAdicLocalField p) {n : ℕ} (hn : ¬ p ∣ n)
    (hn0 : n ≠ 0) (hdvd : n ∣ Nat.card 𝓀[K.carrier] - 1) :
    contHomCard K.absGal n = n * n := by
  obtain ⟨ζ, hζ⟩ := exists_isPrimitiveRoot_mem_carrier K hdvd
  -- `K.absGal` は `AlgebraicClosure K.carrier ≃ₐ[K.carrier] AlgebraicClosure K.carrier` の
  -- `abbrev` なので、Kummer 双対はそのまま当たる(移送は要らない)。
  have h1 : contHomCard K.absGal n
      = ((powMonoidHom n : (K.carrier)ˣ →* (K.carrier)ˣ).range).index :=
    card_contHom_eq_index_powRange K.carrier hn0 hζ
  rw [h1, index_powRange_carrierUnits K hn, Nat.gcd_eq_left hdvd]

end ABC3.Found.PGC
