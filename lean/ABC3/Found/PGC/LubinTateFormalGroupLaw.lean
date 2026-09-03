import ABC3.Found.PGC.LubinTateLimit

/-!
# Lubin-Tate 形式群法則 `F_f`(`g = f` の特殊化、`sorry` 無し)

`Found/PGC/LubinTateLimit.lean::LubinTateF` は一般の `f,g`(両方とも
`≡πX(mod deg2)`)から出発したが、**古典的な Lubin-Tate 形式群法則**
`F_f` はその特殊化 `g := f` にあたる——`f` 自身が `F_f` の**自己準同型**
(乗法 `[π]`)を与える、という古典的な設定に一致する。

## 数学的な内容

`f ∈ 𝔉_π`(`f≡πX(mod deg2)`・剰余体上 `f̄(T)=T^q`)を1つ固定すると、
`F_f(X,Y) := LubinTateF(f,f)` は関数等式

```
F_f(f(X),f(Y)) = f(F_f(X,Y))
```

を満たす——これは「`f` が `F_f` の自己準同型である」という古典的な事実
そのもの(Lubin–Tate 1965、"formal complex multiplication"の出発点)。
`f` はこの意味で乗法元 `π ∈ 𝒪_K` に対応する自己準同型 `[π]_{F_f}` を
与える——`𝒪_K → End(F_f)` という環準同型(formal 𝒪_K-module 構造)の
第一歩。

## まだ無いもの

`F_f` が実際に**形式群法則**であること(単位元則 `F_f(X,0)=X`・
結合律・可換律)は、Lubin-Tate の「一意性補題」(2つの候補が同じ次数1の
係数を持ち同じ関数等式を満たせば等しい、という事実)を経由するのが
標準的だが、本ファイルでは未着手——`coeff_subst_linearize`/
`coeff_subst_g_linearize` と同型の次数ごとの帰納法で示せると見ている。
さらにその先には、`[π^n]_{F_f}` の非零な torsion 点(`f` の反復合成の
根)を使って Lubin-Tate 拡大 `K_{π,n}` を構成し、その上の Galois 作用が
相互律写像そのものを与える、という段が続く——ここまでで局所類体論の
相互律の**構成そのもの**にはまだ距離がある。
-/

namespace ABC3.Found.PGC

variable {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hfres : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff))

include hf0 in
omit [IsLocalRing A] [IsDomain A] [ExpChar (IsLocalRing.ResidueField A) pp]
  [Fintype (IsLocalRing.ResidueField A)] in
/-- `f` を「`g` の役」でも使うための、定数項0の言い換え
(`PowerSeries.coeff 0 f = 0` ↔ `PowerSeries.constantCoeff f = 0`)。 -/
theorem hf0' : PowerSeries.constantCoeff f = 0 := by
  rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0

include hq hπmax hf0 hf1 hfres in
/-- ★★**古典的な Lubin-Tate 形式群法則** `F_f`——`LubinTateF` を `g:=f` で
特殊化したもの。 -/
noncomputable def formalGroupLaw : MvPowerSeries (Fin 2) A :=
  LubinTateF hq hπmax f (hf0' f hf0) hf1 hfres f hf0 hf1 hfres

include hq hπmax hf0 hf1 hfres in
/-- ★★★★★★★★★**`f` は `F_f` の自己準同型である**(Lubin-Tate 形式群法則の
基礎的な性質)。`F_f(f(X),f(Y)) = f(F_f(X,Y))`——`LubinTateF_functional_
equation` を `g:=f` で特殊化するだけで直ちに出る。 -/
theorem formalGroupLaw_f_isEndomorphism :
    MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f)
        (formalGroupLaw hq hπmax f hf0 hf1 hfres) =
      PowerSeries.subst (formalGroupLaw hq hπmax f hf0 hf1 hfres) f :=
  LubinTateF_functional_equation hq hπmax f (hf0' f hf0) hf1 hfres f hf0 hf1 hfres

end ABC3.Found.PGC
