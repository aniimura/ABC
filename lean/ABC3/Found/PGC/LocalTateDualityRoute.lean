import ABC3.Found.PGC.ArtinEquivariance

/-!
# [pGC] Proposition 1.1 —— **2 本目の道**: 局所 Tate 双対性

原文 (pGC p.3, Proposition 1.1):

> The cyclotomic character χ : ΓK → Zp× can be recovered entirely
> group-theoretically from ΓK.

## ★本ファイルの立ち位置

`Found/PGC/ArtinEquivariance.lean` は Proposition 1.1 を **Artin 写像の同変性**
(`ArtinEquivariance`)に帰着させた。本ファイルは**別の道**——**原典が実際に通っている道**
——を敷く。

原典 (pGC p.3) の Proposition 1.1 直前の段落が論拠である:

> Hi(K, M) ≅ H2−i(K, M ∨(1))∨ for i ≥ 0. Here, the "(1)" is a Tate twist,
> and the superscripted "∨"'s denote the "Pontrjagin dual" (i.e., Hom(−, Qp/Zp)).

> if and only if H2(K, M) ≅ Z/pnZ. This is clearly a group-theoretic condition
> on M. Thus, we conclude that the isomorphism class of the ΓK-module Zp(1)
> can be recovered group-theoretically from ΓK.

★★**原典は Artin 同変性を一度も使っていない。** Proposition 1.1 の直後にある局所類体論の
段落 (`ΓKab ≅ (K×)∧`) は **Proposition 1.2 用**であって Proposition 1.1 用ではない。
これは `.txt` 88–135 行の直読で確定した。

## ★本ファイルが**無条件に**証明したこと(仮定ゼロ)

1. **`forall_of_natCard_subtype_eq` / `natCard_subtype_eq_of_forall`**(抽象核 0)——
   有限型の部分型が全体と同じ濃度なら述語は恒真。分岐・付値・Galois・コホモロジーの
   語彙ゼロ。
2. **`natCard_charTwistFixed_eq_iff`**(抽象核 1)——
   2 つの指標 `ψ, χ : G → Rˣ` について、ねじれ `ψ⁻¹χ` の**不動点が全体と同じ濃度**
   ⟺ `ψ = χ`。★これが原典の「`H2(K, M) ≅ Z/pnZ` ⟺ `M ≅ Z/pnZ(1)`」の**代数の中身**
   である(双対性が `|H²|` を `|(M∨(1))^{Γ_K}|` に置き換えた後に残る部分)。
   ★★分岐・付値・Galois・コホモロジーの語彙が **1 語も出ない**。
3. **`cycloCharHom`** —— 円分指標の `mod p^n` 還元は `K.absGal →* (ZMod (p^n))ˣ` という
   **束ねられた指標**であり、その**核は開**(`isOpen_ker_cycloCharHom`)。
   ★これは「`M_ψ` が離散 `Γ_K`-加群である」という双対性の適用条件そのものである。
4. **★★★`cyclotomicCharacterObject_recoverable_of_modPowTransport`** ——
   Proposition 1.1(`ℤ_[p]ˣ` での等式)は **`mod p^n` の合同**に完全に帰着する。
   ★これが本ファイルの主結果である。局所 Tate 双対性は `mod p^n` でしか物を言わないので、
   **目標を双対性が届く高さまで下ろす**この一段が要る。
   逆も証明した(`modPowTransport_of_recoverable`)ので**同値**であり、
   ★**言い換えによって主張を強めていない**ことが確定している。
5. **`LocalTateDualityData.recoverable`** —— 原典の道筋を、`H²` の**位数**が満たすべき
   2 条件(群論性 + 双対性)に切り出し、そこから Proposition 1.1 が出ることを証明した。
   ★`H²` そのものを定義せずに済ませてあるので、mathlib に連続コホモロジーが入った時点で
   `LocalTateDualityData` を**組み立てるだけ**になる。

## ★★仮定に置いたもの(名指し)

**`LocalTateDualityData`(構造体)は仮定である。証明していない。**
その 3 つのフィールドはそれぞれ次を仮定している:

* `cardH2` —— 開核をもつ指標 `ψ` に対する `|H²(Γ_K, M_ψ)|`(の値)の**存在**。
* `isGroupTheoretic` —— `α : Γ_K ≃ₜ* Γ_{K'}` に沿って `|H²|` が**不変**であること
  (原典の「This is clearly a group-theoretic condition on M」)。
* `cardH2_eq_natCard` —— **局所 Tate 双対性**(位数のレベル):
  `|H²(Γ_K, M_ψ)| = |(M_ψ∨(1))^{Γ_K}|`。これは原典が [2] Proposition 3.8 に投げた
  `Hi(K, M) ≅ H2−i(K, M ∨(1))∨` を `i = 2` で使い、Pontrjagin 双対が有限アーベル群の
  位数を保つことを合わせたものである。

★**この構造体を作ることは本ファイルではしていない。** できない理由は下の測定を見よ。

## ★★★mathlib 測定(2026-09-07 実施。★2026-09-04 の記録から**変わった**)

`.cache/mathlib-index.txt`(2026-09-07 09:38 生成)への grep と、
`tools/leanfile.mjs` による型検査で測った。

| 測ったもの | 2026-09-04 の記録 | ★2026-09-07 の実測 |
|---|---|---|
| 離散 `groupCohomology` | 在る(`[Group G]` のみ) | 在る。★`LowDegree`(145 宣言)/ `Functoriality`(63)/ ★**`LongExactSequence`(19)** / ★**`FiniteCyclic`(9)** / ★**`Hilbert90`(7)** / `Shapiro`(2) |
| `TateCohomology` | 在る | 在る(`RepresentationTheory/Homological/TateCohomology/Basic.lean`, 33 宣言) |
| `continuousCohomology` | `RepresentationTheory/...`、対象は `Action (TopModuleCat R) G` | ★**移動していた**。`Algebra/Category/ContinuousCohomology/Basic.lean`。★新たに `RepresentationTheory/Continuous/Basic.lean`(`ContRepresentation`)と `RepresentationTheory/Continuous/TopRep.lean`(`TopRep`)が在る |
| `H²` の inflation-restriction | 測っていない | ★**無い**。`groupCohomology.H1InfRes` / `H1InfRes_exact` は**次数 1 だけ**。`H2InfRes` は `Unknown constant` |
| 有限群の塔から副有限群への colimit | 「長完全列が無い」 | ★**無い**。`groupCohomology.colimitIso` は `Unknown constant`。inflation は `infNatTrans` として自然変換の形で在るが、極限を取る道具は無い |
| `M_ψ = ZMod m(ψ)` を作れるか | 測っていない | ★**作れる**。`Rep.of ((Units.coeHom _).comp ((Units.map (Algebra.lmul _ _).toRingHom.toMonoidHom).comp ψ))` が通る |
| `Nat.card (groupCohomology.H2 A)` | 測っていない | ★**型は付く**(`ModuleCat` の `CoeSort` 経由) |
| `Finite (groupCohomology.H2 A)`(`[Finite G]` 付き) | 測っていない | ★**無い**(`failed to synthesize`)。★**有限群・有限係数ですら `H²` の有限性インスタンスが無い**ので、`|H²| = p^n` を述べるだけでも有限性を手で作る必要がある |
| 群同型に沿った `H²` の同型 | 測っていない | ★**無い**。`e : G ≃* H` について `H2 (Rep.res e.toMonoidHom A) ≅ H2 A` は `exact?` が閉じられない。`groupCohomology.map` と `congr` から組めるはずだが**束ねられていない** |
| Brauer 群 | 不在 | `Algebra/BrauerGroup/Defs.lean` に `BrauerGroup` / `CSA` / `IsBrauerEquivalent` の**定義だけ**(9 宣言)。★**不変写像 `Br(K) ≅ ℚ/ℤ` も crossed product も無い** |
| 局所 Tate 双対性・Poitou–Tate | 不在 | ★**不在のまま**(`tateduality` / `localduality` / `poitou` は 0 件) |
| 局所類体論の相互写像 | 不在 | ★**不在のまま**(`artinReciprocity` / `localReciprocity` / `normResidue` / `invariantMap` は 0 件) |

★★**結論(測定)**: 「有限次 `L/K` の `Gal(L/K)` に離散 `groupCohomology` を直接当てる」道は
**`M_ψ` を作って `Nat.card (H² …)` と書くところまでは行ける**。そこで止まる理由は 3 つあり、
どれも**離散側の欠落**である:
(a) `H²` の**有限性インスタンスが無い**、
(b) **群同型に沿った `H²` の同型が束ねられていない**、
(c) **`H²` の inflation-restriction と塔の colimit が無い**ので、
有限次で得た情報を `Γ_K` に上げられない。
★(c) が本質的な壁である。(a)(b) は分量の問題にすぎない。

## ★残りを「存在量化子を含まない 1 文」に縮めたもの

**`CyclotomicCharacterModPowTransport p`**:

> 任意の `p` 進局所体 `K`, `K'`、任意の位相群同型 `α : Γ_K ≃ₜ* Γ_{K'}`、
> 任意の `n : ℕ`、任意の `g : Γ_K` について
> `χ_{K'} (α g) mod p^n = χ_K g mod p^n`。

★**存在量化子を含まない**(既存の接続点 `cyclotomicCharacterObject_transport_of_moduleEquiv`
の仮説は `∃ φ : ℤ_[p] ≃ₗ[ℤ_[p]] ℤ_[p]` を含んでいた)。
★`recoverable_iff_modPowTransport` により **Proposition 1.1 と同値**である。

## 逸脱の記録

* **逸脱 1**: 原典は `H2(K, M) ≅ Z/pnZ`(**同型**)と書く。本ファイルは
  `|H²| = p^n`(**位数の一致**)で置き換えた。`M` は位数 `p^n` の巡回群なので
  「位数 `p^n` の部分群 ⟺ 全体」の形で使うぶんには等価であり、後続に影響しない。
  ★位数に落としたのは、`H²` を定義せずに「双対性が言うこと」だけを取り出すためである。
* **逸脱 2**: 原典の「continuous ΓK-action」を、本ファイルは
  「指標であって核が開」(`ψ : Γ_K →* (ZMod (p^n))ˣ` かつ `IsOpen ψ.ker`)で表した。
  `(ZMod (p^n))ˣ` に位相インスタンスが無い(実測: `failed to synthesize
  TopologicalSpace (ZMod (p ^ n))ˣ`)ので `Continuous` が書けないためである。
  離散有限群への準同型については両者は同値であり、★**仮定を弱める方向**の書き換えである。
* **逸脱 3**: `M ∨(1)` の Pontrjagin 双対と Tate 捩れを、指標のレベルで
  `ψ⁻¹ · χ` と書いた(`charTwistFixed`)。階数 1 なのでこれで尽きている。

## 退化の自己検査

* ★`natCard_charTwistFixed_self` —— `ψ = χ` のとき不動点はちょうど `p^n` 個。
  **壁は空虚ではない。**
* ★`recoverable_iff_modPowTransport` —— 残った 1 文は Proposition 1.1 と**同値**。
  弱めすぎても強めすぎてもいない。
* ★`LocalTateDualityData` は「双対性がある」と仮定して配線だけ書いたものでは**ない**。
  `cardH2_eq_natCard` は**位数の等式**という検証可能な形であり、そこから
  `cardH2_eq_iff`(判定条件)は**定理として導いてある**(抽象核 1 による)。
  すなわち「結論を仮定に書いた」形にはなっていない。
* ★Kummer 理論は使っていない(`µ_n` を係数に要求するので円分子の構成には使えない)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC

universe u

/-! ## 1. 抽象核 —— 分岐・付値・Galois・コホモロジーの語彙ゼロ

★この節には `PAdicLocalField` も `cyclotomicCharacter` も出てこない。
一般の有限型と一般の指標の話である。 -/

section AbstractCore

variable {M : Type*} {P : M → Prop}

/-- **★抽象核 0(前半)** —— 有限型 `M` の部分型 `{x // P x}` が `M` と同じ濃度なら、
`P` は恒真。

★分岐・付値・Galois・コホモロジーの語彙ゼロ。使うのは
「単射 + 濃度一致 ⟹ 全射」(`Nat.bijective_iff_injective_and_card`)だけである。

★この一行が原典の「`H2(K, M) ≅ Z/pnZ` ⟹ `M ≅ Z/pnZ(1)`」の向きを担う:
双対性が `|H²|` を「ねじれの不動点の個数」に翻訳したあと、
「不動点が全部 ⟹ 作用が自明」を出すのがこれである。 -/
theorem forall_of_natCard_subtype_eq [Finite M]
    (h : Nat.card {x : M // P x} = Nat.card M) : ∀ x, P x := by
  have hbij : Function.Bijective (Subtype.val : {x : M // P x} → M) :=
    (Nat.bijective_iff_injective_and_card _).2 ⟨Subtype.val_injective, h⟩
  intro x
  obtain ⟨⟨y, hy⟩, rfl⟩ := hbij.2 x
  exact hy

/-- **★抽象核 0(後半)** —— 逆向き。`P` が恒真なら部分型は全体と同じ濃度。

★`Finite M` は不要である(`Nat.card` は無限でも定義されている)。 -/
theorem natCard_subtype_eq_of_forall (h : ∀ x, P x) :
    Nat.card {x : M // P x} = Nat.card M :=
  Nat.card_eq_of_bijective Subtype.val ⟨Subtype.val_injective, fun x => ⟨⟨x, h x⟩, rfl⟩⟩

variable {G R : Type*} [Monoid R] {ψ χ : G → Rˣ}

/-- **指標のねじれの不動点**。

`ψ, χ : G → Rˣ` に対し、`R` を `g • x = (ψ g)⁻¹ (χ g) x` で `G`-加群と見たときの不動点。

★原典の記号では `M_ψ∨(1)` の `Γ_K`-不動点、すなわち `H0(K, M∨(1))` にあたる。
階数 1 なので「`ψ` を `χ` でねじる」ことがそのまま Tate 捩れと Pontrjagin 双対の合成になる。 -/
def charTwistFixed (ψ χ : G → Rˣ) : Set R :=
  {x : R | ∀ g : G, (((ψ g)⁻¹ * χ g : Rˣ) : R) * x = x}

/-- **★★抽象核 1** —— ねじれの不動点が全体と同じ濃度 ⟺ 2 つの指標が一致。

★★これが原典の

> if and only if H2(K, M) ≅ Z/pnZ

の**代数の中身**である。局所 Tate 双対性が `|H²(K, M_ψ)|` を
`|(M_ψ∨(1))^{Γ_K}|` に置き換えたあと、残るのはこの初等的な同値だけになる。

★分岐・付値・Galois・コホモロジーの語彙が **1 語も出ない**。
`G` は単なる型でよく(群である必要すらない)、`R` は有限モノイドでよい。 -/
theorem natCard_charTwistFixed_eq_iff [Finite R] :
    Nat.card (charTwistFixed ψ χ) = Nat.card R ↔ ψ = χ := by
  constructor
  · intro h
    have h1 := forall_of_natCard_subtype_eq h 1
    funext g
    have hone : ((ψ g)⁻¹ * χ g : Rˣ) = 1 := Units.ext (by simpa using h1 g)
    rwa [inv_mul_eq_one] at hone
  · rintro rfl
    exact natCard_subtype_eq_of_forall (fun x g => by simp)

/-- **★抽象核 1'** —— 濃度を外から与える形。`Nat.card R` を `p ^ n` に読み替えるときに使う。

★`rw [← hm]` が動くよう `m` を明示引数で持つ(`lean-idioms.md` #150)。 -/
theorem natCard_charTwistFixed_eq_iff' [Finite R] {m : ℕ} (hm : Nat.card R = m) :
    Nat.card (charTwistFixed ψ χ) = m ↔ ψ = χ := by
  rw [← hm]
  exact natCard_charTwistFixed_eq_iff

end AbstractCore

/-! ## 2. 具体層 —— 円分指標の `mod p^n` 還元

★ここから `PAdicLocalField` が出てくる。抽象核とは別の宣言にしてある。 -/

section ModPow

variable {p : ℕ} [Fact p.Prime]

/-- **円分指標の `mod p^n` 還元**(値だけ)。

`χ_K : Γ_K → ℤ_[p]ˣ` を `ℤ_[p] → ZMod (p^n)` で押し出したもの。 -/
noncomputable def cycloCharUnitsModPow (K : PAdicLocalField p) (n : ℕ) (g : K.absGal) :
    (ZMod (p ^ n))ˣ :=
  Units.map (PadicInt.toZModPow (p := p) n : ℤ_[p] →+* ZMod (p ^ n)).toMonoidHom
    (cyclotomicCharacter K.closure p g.toRingEquiv)

/-- `cycloCharUnitsModPow` の台は `χ_K` の `toZModPow` そのもの。 -/
theorem coe_cycloCharUnitsModPow (K : PAdicLocalField p) (n : ℕ) (g : K.absGal) :
    ((cycloCharUnitsModPow K n g : (ZMod (p ^ n))ˣ) : ZMod (p ^ n))
      = PadicInt.toZModPow n (cyclotomicCharacter K.closure p g.toRingEquiv : ℤ_[p]) :=
  rfl

/-- **円分指標の `mod p^n` 還元**(束ねた準同型)。

★`cyclotomicCharacter K.closure p : (K.closure ≃+* K.closure) →* ℤ_[p]ˣ` は mathlib で
既に `MonoidHom` であり、`MulSemiringAction.toRingAut K.absGal K.closure` が
`Γ_K = Gal(K̄/K)` から `RingAut` への準同型を与える。両者の合成に `Units.map` を被せる。 -/
noncomputable def cycloCharHom (K : PAdicLocalField p) (n : ℕ) :
    K.absGal →* (ZMod (p ^ n))ˣ :=
  (Units.map (PadicInt.toZModPow (p := p) n : ℤ_[p] →+* ZMod (p ^ n)).toMonoidHom).comp
    ((cyclotomicCharacter K.closure p).comp (MulSemiringAction.toRingAut K.absGal K.closure))

@[simp]
theorem cycloCharHom_apply (K : PAdicLocalField p) (n : ℕ) (g : K.absGal) :
    cycloCharHom K n g = cycloCharUnitsModPow K n g :=
  rfl

theorem coe_cycloCharHom (K : PAdicLocalField p) (n : ℕ) :
    ⇑(cycloCharHom K n) = cycloCharUnitsModPow K n :=
  rfl

/-- `muFixer K (p^n)`(`p^n` 乗根をすべて止める元)は `χ mod p^n` の核に入る。

★`Found/PGC/CyclotomicRecoverable.lean::toZModPow_cyclotomicCharacter_eq_one_of_mem_muFixer`
の言い換え。 -/
theorem muFixer_le_ker_cycloCharHom (K : PAdicLocalField p) (n : ℕ) :
    muFixer K (p ^ n) ≤ (cycloCharHom K n).ker := by
  intro g hg
  simp only [MonoidHom.mem_ker]
  refine Units.ext ?_
  rw [cycloCharHom_apply, coe_cycloCharUnitsModPow, Units.val_one]
  exact toZModPow_cyclotomicCharacter_eq_one_of_mem_muFixer K n hg

/-- **★`χ mod p^n` の核は開**。

★これが「`M_ψ` が離散 `Γ_K`-加群である」——原典の "continuous ΓK-action"——の中身である。
`muFixer K (p^n)` が開である(`isOpen_muFixer`)ことと、
**開部分群を含む部分群は開**(`Subgroup.isOpen_mono`)から出る。

★`(ZMod (p^n))ˣ` には位相インスタンスが無い(実測)ので `Continuous` とは書けない。
離散有限群への準同型については「核が開」と「連続」は同値なので、こちらで表す。 -/
theorem isOpen_ker_cycloCharHom (K : PAdicLocalField p) (n : ℕ) :
    IsOpen (((cycloCharHom K n).ker : Subgroup K.absGal) : Set K.absGal) :=
  Subgroup.isOpen_mono (muFixer_le_ker_cycloCharHom K n) (isOpen_muFixer K (p ^ n))

/-- 開核をもつ指標を位相群同型で引き戻しても核は開。

★`α` の連続性しか使わない配管。 -/
theorem isOpen_ker_comp (K K' : PAdicLocalField p)
    (α : ContinuousMulEquiv K.absGal K'.absGal) (n : ℕ)
    (ψ : K'.absGal →* (ZMod (p ^ n))ˣ)
    (h : IsOpen ((ψ.ker : Subgroup K'.absGal) : Set K'.absGal)) :
    IsOpen (((ψ.comp (α : K.absGal ≃* K'.absGal).toMonoidHom).ker : Subgroup K.absGal)
      : Set K.absGal) := by
  rw [← MonoidHom.comap_ker, Subgroup.coe_comap]
  exact h.preimage α.continuous

end ModPow

/-! ## 3. ★★★残った 1 文(存在量化子なし)と、Proposition 1.1 との同値 -/

section ModPowTransport

variable {p : ℕ} [Fact p.Prime]

def CyclotomicCharacterModPowTransport.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- **★★★局所 Tate 双対性の道に残った 1 文**(存在量化子を含まない)。

原文 (pGC p.3):
> The cyclotomic character χ : ΓK → Zp× can be recovered entirely
> group-theoretically from ΓK.

を `mod p^n` に下ろしたもの:

> 任意の `K`, `K'`, `α : Γ_K ≃ₜ* Γ_{K'}`, `n`, `g` について
> `χ_{K'} (α g) ≡ χ_K g  (mod p^n)`。

★★**なぜ下ろすのか**: 局所 Tate 双対性は「有限長の `ℤ_p`-加群 `M`」についての主張なので、
`ℤ_[p]ˣ` での等式には直接届かない。`mod p^n` まで下ろすと、原典の段落がそのまま使える
高さになる。

★★`recoverable_iff_modPowTransport` により **Proposition 1.1 と同値**である
(弱めても強めてもいない)。 -/
def CyclotomicCharacterModPowTransport (p : ℕ) [Fact p.Prime] : Prop :=
  ∀ (K K' : PAdicLocalField p) (α : ContinuousMulEquiv K.absGal K'.absGal)
    (n : ℕ) (g : K.absGal),
    cycloCharUnitsModPow K' n (α g) = cycloCharUnitsModPow K n g

def cyclotomicCharacterObject_recoverable_of_modPowTransport.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- **★★★本ファイルの主結果** —— `mod p^n` の合同がすべての `n` で成り立てば
Proposition 1.1 が出る。

★`ℤ_[p]` が射影極限であること(`PadicInt.ext_of_toZModPow`)だけを使う。
★**仮定ゼロ**。局所 Tate 双対性も Artin 同変性も使っていない。 -/
theorem cyclotomicCharacterObject_recoverable_of_modPowTransport
    (h : CyclotomicCharacterModPowTransport p) :
    (cyclotomicCharacterObject (p := p)).RecoverableFromAbsGal := by
  refine cyclotomicCharacterObject_recoverable_iff.mpr ?_
  intro K K' α g
  refine Units.ext (PadicInt.ext_of_toZModPow.mp (fun n => ?_))
  have hn := congrArg Units.val (h K K' α n g)
  rwa [coe_cycloCharUnitsModPow, coe_cycloCharUnitsModPow] at hn

/-- **逆向き** —— Proposition 1.1 から `mod p^n` の合同が出る(自明)。

★これがあるので「`mod p^n` に下ろす」書き換えは**主張を強めていない**。
退化の自己検査である。 -/
theorem modPowTransport_of_recoverable
    (h : (cyclotomicCharacterObject (p := p)).RecoverableFromAbsGal) :
    CyclotomicCharacterModPowTransport p := by
  intro K K' α n g
  have hg := cyclotomicCharacterObject_recoverable_iff.mp h K K' α g
  unfold cycloCharUnitsModPow
  rw [hg]

/-- **★★残った 1 文は Proposition 1.1 と同値**。 -/
theorem recoverable_iff_modPowTransport :
    (cyclotomicCharacterObject (p := p)).RecoverableFromAbsGal
      ↔ CyclotomicCharacterModPowTransport p :=
  ⟨modPowTransport_of_recoverable, cyclotomicCharacterObject_recoverable_of_modPowTransport⟩

/-- `mod p^n` の合同から `ℤ_[p]ˣ` での等式。 -/
theorem cyclotomicCharacter_eq_of_modPowTransport
    (h : CyclotomicCharacterModPowTransport p) (K K' : PAdicLocalField p)
    (α : ContinuousMulEquiv K.absGal K'.absGal) (g : K.absGal) :
    cyclotomicCharacter K'.closure p (α g).toRingEquiv
      = cyclotomicCharacter K.closure p g.toRingEquiv :=
  cyclotomicCharacterObject_recoverable_iff.mp
    (cyclotomicCharacterObject_recoverable_of_modPowTransport h) K K' α g

/-- **★既存の接続点への橋**。

`Found/PGC/CyclotomicRecovery.lean::cyclotomicCharacterObject_transport_of_moduleEquiv` は

  `∃ φ : ℤ_[p] ≃ₗ[ℤ_[p]] ℤ_[p], ∀ g x, φ (χ_K g • x) = χ_{K'} (α g) • φ x`

を仮説に取る。★**双対性の道はこの仮説を直接は作らない**(双対性が言うのは有限長加群
`M ≅ ℤ/p^n` についてであって、`ℤ_[p]` 上の線形同型ではない)。
本補題は、`mod p^n` の合同からならこの仮説が**作れる**ことを示す(`φ = id` でよい)。

★すなわち `CyclotomicCharacterModPowTransport` は既存の接続点の仮説より**強くない**——
むしろ `∃` を消したぶん扱いやすい。 -/
theorem moduleEquiv_hypothesis_of_modPowTransport
    (h : CyclotomicCharacterModPowTransport p) :
    ∀ (K K' : PAdicLocalField p) (α : ContinuousMulEquiv K.absGal K'.absGal),
      ∃ φ : ℤ_[p] ≃ₗ[ℤ_[p]] ℤ_[p], ∀ (g : K.absGal) (x : ℤ_[p]),
        φ ((cyclotomicCharacter K.closure p g.toRingEquiv : ℤ_[p]) • x)
          = (cyclotomicCharacter K'.closure p (α g).toRingEquiv : ℤ_[p]) • φ x := by
  intro K K' α
  refine ⟨LinearEquiv.refl _ _, fun g x => ?_⟩
  rw [cyclotomicCharacter_eq_of_modPowTransport h K K' α g]
  rfl

/-- **★2 本の道は同じ 1 点に着地する**。

`Found/PGC/ArtinEquivariance.lean` の `ArtinEquivariance`(Lubin–Tate 経由の道)からも
`CyclotomicCharacterModPowTransport` が出る。★どちらか一方が埋まれば Proposition 1.1 が閉じる。 -/
theorem modPowTransport_of_artinEquivariance (h : ArtinEquivariance p) :
    CyclotomicCharacterModPowTransport p :=
  modPowTransport_of_recoverable (cyclotomicCharacter_recoverable_of_artinEquivariance h)

end ModPowTransport

/-! ## 4. ★★原典の道筋 —— 局所 Tate 双対性を「位数の 2 条件」に切り出す

★★**この節の `LocalTateDualityData` は仮定である。証明していない。**
mathlib に副有限群の連続 `H²` も局所 Tate 双対性も無い(上の測定表を見よ)。

★ただし**結論を仮定に書いた形にはしていない**。フィールドは
「群論性」と「位数の等式」という 2 つの検証可能な条件だけであり、
そこから判定条件 `cardH2_eq_iff` は**定理として**導いてある(抽象核 1 による)。 -/

section Duality

variable {p : ℕ} [Fact p.Prime]

def LocalTateDualityData.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- **★★局所 Tate 双対性が Proposition 1.1 に供給すべきもの**(仮定の束)。

原文 (pGC p.3):
> Hi(K, M) ≅ H2−i(K, M ∨(1))∨ for i ≥ 0. Here, the "(1)" is a Tate twist,
> and the superscripted "∨"'s denote the "Pontrjagin dual" (i.e., Hom(−, Qp/Zp)).

> if and only if H2(K, M) ≅ Z/pnZ. This is clearly a group-theoretic condition
> on M. Thus, we conclude that the isomorphism class of the ΓK-module Zp(1)
> can be recovered group-theoretically from ΓK.

フィールドの意味:

* `cardH2 K n ψ` —— `|H²(Γ_K, M_ψ)|`。ここで `M_ψ` は `ZMod (p^n)` に `Γ_K` が指標 `ψ` で
  作用する加群(原典の `M`)。★`H²` そのものは定義していない。**位数だけ**を持つ。
* `isGroupTheoretic` —— 原典の「This is clearly a group-theoretic condition on M」。
  `α : Γ_K ≃ₜ* Γ_{K'}` に沿って `|H²|` が不変であること。
* `cardH2_eq_natCard` —— **局所 Tate 双対性**を位数のレベルで書いたもの。
  原典の `i = 2` の場合 `H2(K, M) ≅ H0(K, M∨(1))∨` と、
  有限アーベル群の Pontrjagin 双対が位数を保つことを合わせると
  `|H²(Γ_K, M_ψ)| = |(M_ψ∨(1))^{Γ_K}|` になる。右辺が `charTwistFixed` である。

★★**仮定を置いたのはこの 3 つだけである。** これを実際に作るには mathlib に
(a) 副有限群の連続 `H²`、(b) 局所 Tate 双対性、(c) `H²` の有限性、が要る。
2026-09-07 の実測ではどれも無い。 -/
structure LocalTateDualityData (p : ℕ) [Fact p.Prime] where
  /-- `|H²(Γ_K, M_ψ)|`。`ψ` は `mod p^n` の指標。 -/
  cardH2 : ∀ (K : PAdicLocalField p) (n : ℕ), (K.absGal →* (ZMod (p ^ n))ˣ) → ℕ
  /-- 原典の「This is clearly a group-theoretic condition on M」。 -/
  isGroupTheoretic : ∀ (K K' : PAdicLocalField p)
    (α : ContinuousMulEquiv K.absGal K'.absGal) (n : ℕ)
    (ψ : K'.absGal →* (ZMod (p ^ n))ˣ),
    IsOpen ((ψ.ker : Subgroup K'.absGal) : Set K'.absGal) →
    cardH2 K n (ψ.comp (α : K.absGal ≃* K'.absGal).toMonoidHom) = cardH2 K' n ψ
  /-- **局所 Tate 双対性**(位数のレベル): `|H²(Γ_K, M_ψ)| = |(M_ψ∨(1))^{Γ_K}|`。 -/
  cardH2_eq_natCard : ∀ (K : PAdicLocalField p) (n : ℕ)
    (ψ : K.absGal →* (ZMod (p ^ n))ˣ),
    IsOpen ((ψ.ker : Subgroup K.absGal) : Set K.absGal) →
    cardH2 K n ψ = Nat.card (charTwistFixed (⇑ψ) (cycloCharUnitsModPow K n))

namespace LocalTateDualityData

/-- **★★原典の判定条件は定理である**(仮定ではない)。

> M is isomorphic as a ΓK-module to Z/pnZ(1) if and only if H2(K, M) ≅ Z/pnZ

★双対性(`cardH2_eq_natCard`)と抽象核 1(`natCard_charTwistFixed_eq_iff'`)だけから出る。
★**ここが「結論を仮定に書いていない」ことの証拠**である。 -/
theorem cardH2_eq_iff (D : LocalTateDualityData p) (K : PAdicLocalField p) (n : ℕ)
    (ψ : K.absGal →* (ZMod (p ^ n))ˣ)
    (hψ : IsOpen ((ψ.ker : Subgroup K.absGal) : Set K.absGal)) :
    D.cardH2 K n ψ = p ^ n ↔ ⇑ψ = cycloCharUnitsModPow K n := by
  rw [D.cardH2_eq_natCard K n ψ hψ]
  exact natCard_charTwistFixed_eq_iff' (Nat.card_zmod _)

/-- **★★原典の論証をそのまま辿る一段**。

`α : Γ_K ≃ₜ* Γ_{K'}` を取る。`χ_{K'}` は `K'` 側で判定条件を満たす(`cardH2 = p^n`)。
`|H²|` は群論的なので、`α` で引き戻した `χ_{K'} ∘ α` も `K` 側で `cardH2 = p^n` を満たす。
再び判定条件を `K` 側で使うと `χ_{K'} ∘ α = χ_K`。 -/
theorem modPowTransport (D : LocalTateDualityData p) :
    CyclotomicCharacterModPowTransport p := by
  intro K K' α n g
  have hopen : IsOpen (((cycloCharHom K' n).ker : Subgroup K'.absGal) : Set K'.absGal) :=
    isOpen_ker_cycloCharHom K' n
  have h1 : D.cardH2 K' n (cycloCharHom K' n) = p ^ n :=
    (D.cardH2_eq_iff K' n _ hopen).2 (coe_cycloCharHom K' n)
  have h2 : D.cardH2 K n
      ((cycloCharHom K' n).comp (α : K.absGal ≃* K'.absGal).toMonoidHom) = p ^ n := by
    rw [D.isGroupTheoretic K K' α n _ hopen]
    exact h1
  exact congrFun ((D.cardH2_eq_iff K n _ (isOpen_ker_comp K K' α n _ hopen)).1 h2) g

def recoverable.src : ABC3.Meta.Source :=
  { paper := "pGC", pdfPage := 3, item := "Proposition 1.1", sectionId := "prop-1-1" }

/-- **★★★[pGC] Proposition 1.1 —— 局所 Tate 双対性の道**。

原文 (pGC p.3):
> The cyclotomic character χ : ΓK → Zp× can be recovered entirely
> group-theoretically from ΓK.

★★`LocalTateDualityData` は**仮定**である(本ファイルでは作っていない)。
これが得られれば Proposition 1.1 が閉じる、というのが本定理の内容である。 -/
theorem recoverable (D : LocalTateDualityData p) :
    (cyclotomicCharacterObject (p := p)).RecoverableFromAbsGal :=
  cyclotomicCharacterObject_recoverable_of_modPowTransport D.modPowTransport

end LocalTateDualityData

/-! ### 退化の自己検査 -/

/-- **★壁は空虚ではない** —— `ψ = χ` のとき、ねじれの不動点はちょうど `p^n` 個ある。

★`natCard_cyclotome_eq`(位数ちょうど `p^n`)と同じ役割の検査である。
`charTwistFixed` が空だったり全体だったりして自明になっていないことを確かめている。 -/
theorem natCard_charTwistFixed_self (K : PAdicLocalField p) (n : ℕ) :
    Nat.card (charTwistFixed (cycloCharUnitsModPow K n) (cycloCharUnitsModPow K n)) = p ^ n :=
  (natCard_charTwistFixed_eq_iff' (Nat.card_zmod _)).2 rfl

end Duality

end ABC3.Found.PGC
