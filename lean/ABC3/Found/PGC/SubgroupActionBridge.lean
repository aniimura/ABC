import ABC3.Found.PGC.HasseArf
import ABC3.Found.PGC.HasseArfStrongInduction

/-!
# 部分群 `H` を群として見る橋 —— Hasse-Arf(Theorem 6.11)の段 1+2 と段 3 の接続

★★**本ノードは Yoshida 2008 Theorem 6.11 の「部分群と群の橋」のみを担当する。**
段 1(`G` 巡回)は `Found/PGC/HasseArf.lean`
(`exists_natCast_herbrandPhiGroup_of_isCyclic`)、
段 2(`G = G_1`)は `Found/PGC/HasseArfStrongInduction.lean`
(`exists_natCast_herbrandPhiGroup_of_lowerRamificationGroup_one_eq_top`)、
段 3(`G ≠ G_1`)は `Found/PGC/HasseArf.lean`
(`exists_natCast_herbrandPhiGroup_of_tame_quotient`)にある。
★本ファイルは**既存の 3 ファイルを 1 行も書き換えない**(import のみ)。

## ★なぜ橋が要るか

木には `φ` が 2 つある。

* `herbrandPhi α (H : Subgroup G) n` —— `G` の**中の部分群** `H` に対する `φ_H`。
* `herbrandPhiGroup G α n` —— **群 `G` そのもの**に対する `φ_G`。

段 1・段 2 の結論は `herbrandPhiGroup`(群そのもの)の形であり、
段 3 が要求する仮定 `hk` は `herbrandPhi π G_1`(部分群としての `G_1`)の形である。
両者を繋ぐには `↥H` を**群として** `B` に作用させ、

    herbrandPhi α H n = herbrandPhiGroup ↥H α n

を示す必要がある。それが `herbrandPhi_eq_herbrandPhiGroup_subtype`(§3)である。

## ★★在庫の訂正 —— 4 本のうち 3 本は mathlib に既にあった

先行ノード(Y21)は「`MulSemiringAction ↥H B` は mathlib にインスタンスが無い」と
報告していたが、**これは誤りである**。実測(`.cache/mathlib-index.txt` を名前空間で grep)で
次の 3 本が見つかり、いずれも `inferInstance` で解決する:

| 要るもの | mathlib の宣言 | 場所 |
|---|---|---|
| `MulSemiringAction ↥H B` | `Subgroup.mulSemiringAction` | `Mathlib/Algebra/Ring/Action/Subobjects.lean:40` |
| `SMulCommClass ↥H A B` | `Subgroup.smulCommClass_left` | `Mathlib/Algebra/Group/Subgroup/Actions.lean:43` |
| `FaithfulSMul ↥H B` | 無名 instance | `Mathlib/Algebra/Group/Subgroup/Actions.lean:59` |

さらに `lowerRamificationGroup` の一致も mathlib の

    AddSubgroup.subgroupOf_inertia :
      (I.inertia G).subgroupOf H = I.inertia ↥H       -- Mathlib/Algebra/Group/Subgroup/Basic.lean:1077

がそのまま与える(`Ideal.inertia` は `(Submodule.toAddSubgroup I).inertia` の
reducible な別名なので、`lowerRamificationGroup` の定義を展開する必要すら無い)。
★本ファイルが新たに定義する `instance` は **1 つも無い**。
§1 の `example` 3 本は「本当に解決するか」の回帰検査である。

## ★`⊤` の場合(`herbrandPhiGroup_eq_herbrandPhi_top`)との比較

既にあった `⊤` 版は `Fintype.sum_equiv Subgroup.topEquiv` と `Nat.card_congr` を
経由していた(`↥(⊤ : Subgroup G)` と `G` は別の型なので和を移す必要がある)。
★**一般の `H` では逆に易しく、`rfl` で閉じる。**
`herbrandPhi α H` も `herbrandPhiGroup ↥H α` も同じ添字型 `↥H` の上の `phiOf` であり、
`τ : ↥H` に対する `ramIndex (G := ↥H) α τ` と `ramIndex (G := G) α ↑τ` は
定義的に等しい(`ramIndex_subtype`)からである。

## 退化の自己検査

* ★**`FaithfulSMul ↥H B` を落とすと `ramIndex` が分岐を測らなくなる**
  (`i(σ) = ⊤` が `σ = 1` を意味しなくなる)。§1 の `example` で
  `[FaithfulSMul G B]` から実際に移送できることを確かめてある。
* ★**`[Fintype ↥H]` を落とすと `Nat.card H = 0` で `φ_H` の除算が壊れる**。
  本ファイルの `herbrandPhi` / `herbrandPhiGroup` はすべて `[Fintype ↥H]` を
  明示的に持ち回っており、§4 以降では `[Fintype ↥(lowerRamificationGroup B G 1)]` を
  **インスタンス引数として受け取る**(段 3 の
  `exists_natCast_herbrandPhiGroup_of_tame_quotient` と同じ流儀)。
  ★Y12 が実測したとおり `Fintype ↥(⊤ : Subgroup G)` は `[Fintype G]` から推論されない。
  一般の `H` でも同様なので、推論に頼らず引数にしてある。
* ★**`lowerRamificationGroup B ↥H i = (lowerRamificationGroup B G i).subgroupOf H` は
  `i : ℕ` について全域で成り立つ**(端点でずれない)。理由は上のとおり
  `AddSubgroup.subgroupOf_inertia` がイデアル `𝔪^(i+1)` について一様だからである。
  `i` に関する場合分けも `i = 0` の例外も無い。
* ★`ℕ∞` の切り詰め引き算・除算は 1 つも書いていない(`lean-idioms.md` #102)。
  除算は `phiOf` の定義の中にしか無い。

## ★★`hne` の移送は `n ≥ 1` でしか成り立たない(逸脱ではなく数学)

`G_n ≠ G_{n+1}` から `H_n ≠ H_{n+1}`(`H = G_1`)を出すには
`G_n ≤ H` と `G_{n+1} ≤ H`、すなわち `1 ≤ n` が要る。
`n = 0` では `G_0 = ⊤` と `G_1 = H` の像がどちらも `⊤` に潰れるので移送できない。
★原文もそこは別扱いにしている:

> If n = 0 then φH(0) = 0.

`exists_natCast_herbrandPhi_lowerRamificationGroup_one`(§4)はこの場合分けを
そのまま持っており、`n = 0` は `herbrandPhi_zero` で片付く。

## 逸脱の記録

1. **§6 の `hC` は `G_1 ≠ ⊤` の場合にだけ要求する形にしてある。**
   段 3 が要る `C`(順分岐商 `K′^{G_1}/K` に対応する DVR)と 3 条件
   `hcomp` / `htopC` / `honeC` は、`G = G_1` の場合には不要だからである。
   原典は場合分けの中でだけ `H = G_1` と置いており、本ファイルの形はそれに対応する。
   ★この 3 条件を実際に構成するのは別ノード(順分岐商の塔)であり、
   本ファイルは**それを仮定として受け取ったまま**である。
2. §4・§5 は段 2 が要求する底の設定(`A` が DVR、`IsNoetherian A B`、
   `hresA`、`hAinj`)を引き継ぐ。段 3 単独よりも仮定が増えているが、
   これは段 2 を呼ぶ以上避けられない(原典 §6.1 の設定そのもの)。
-/

namespace ABC3.Found.PGC

open IsLocalRing IsDiscreteValuationRing

def exists_natCast_herbrandPhi_lowerRamificationGroup_one.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 16, item := "Theorem 6.11", sectionId := "thm-6-11" }

def exists_natCast_herbrandPhiGroup_of_abelian.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 16, item := "Theorem 6.11", sectionId := "thm-6-11" }

/-! ## §1 抽象核(純群論)

★分岐・付値・Herbrand 関数の語彙が 1 つも出てこない。 -/

namespace AbstractCore

/-- ★**`Subgroup.subgroupOf` は `H` 以下の部分群の上で単射**。

`K.subgroupOf H = Subgroup.comap H.subtype K` は一般には単射でないが、
`K ≤ H` に制限すれば `Subgroup.map H.subtype` が逆写像になる。
★段 3 に渡す `hne` の移送(`§2`)で使う唯一の群論的事実である。 -/
theorem subgroupOf_injOn {G : Type*} [Group G] {H K L : Subgroup G}
    (hK : K ≤ H) (hL : L ≤ H) (h : K.subgroupOf H = L.subgroupOf H) : K = L := by
  rw [← Subgroup.map_subgroupOf_eq_of_le hK, ← Subgroup.map_subgroupOf_eq_of_le hL, h]

/-- 可換性は部分群に落ちる。★`Subtype.ext` を 1 回使うだけ。 -/
theorem comm_subtype {G : Type*} [Group G] (habel : ∀ x y : G, x * y = y * x) (H : Subgroup G) :
    ∀ x y : ↥H, x * y = y * x := fun x y => Subtype.ext (habel (x : G) (y : G))

end AbstractCore

open AbstractCore

/-! ## §1b インスタンスの移送 —— ★すべて mathlib に既にある

★新しい `instance` は 1 つも宣言しない。以下は「解決すること」の回帰検査である。
`inferInstance` が落ちるようになったら、ここが最初に赤くなる。 -/

/-- `MulSemiringAction ↥H B` —— `Subgroup.mulSemiringAction`(mathlib)。 -/
example {B : Type*} [CommRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    (H : Subgroup G) : MulSemiringAction ↥H B := H.mulSemiringAction

/-- `SMulCommClass ↥H A B` —— `Subgroup.smulCommClass_left`(mathlib)。
★これを落とすと `mem_lowerRamificationGroup_iff`(`B = A[α]` の橋)が使えない。 -/
example {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] {G : Type*} [Group G]
    [MulSemiringAction G B] [SMulCommClass G A B] (H : Subgroup G) :
    SMulCommClass ↥H A B := inferInstance

/-- `FaithfulSMul ↥H B` —— mathlib の無名 instance。
★★**退化検査**: これを落とすと `ramIndex` が分岐を測らなくなる。 -/
example {B : Type*} [CommRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [FaithfulSMul G B] (H : Subgroup G) : FaithfulSMul ↥H B := inferInstance

/-! ## §2 制限の核 —— 下付き分岐群と `i(σ)`

★Herbrand 関数も除算も出てこない。付値が出るのは `ramIndex` の型だけである。 -/

/-- ★★**部分群を群として見たときの下付き分岐群**

    (G_n としての ↥H の分岐群) = (G の分岐群 G_n) の H への制限。

`lowerRamificationGroup B G n = (𝔪_B^(n+1)).inertia G` という定義と
mathlib の `AddSubgroup.subgroupOf_inertia` から**そのまま**出る。

★`n : ℕ` について**全域**で成り立つ(端点の例外が無い)。 -/
theorem lowerRamificationGroup_subtype {B : Type*} [CommRing B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] (H : Subgroup G) (n : ℕ) :
    lowerRamificationGroup B ↥H n = (lowerRamificationGroup B G n).subgroupOf H :=
  (AddSubgroup.subgroupOf_inertia _ H).symm

/-- 所属の形。`σ ∈ H` に対し「`↥H` の `G_n` に入る」と「`G` の `G_n` に入る」は同値。 -/
theorem mem_lowerRamificationGroup_subtype {B : Type*} [CommRing B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {H : Subgroup G} {n : ℕ} {σ : ↥H} :
    σ ∈ lowerRamificationGroup B ↥H n ↔ (σ : G) ∈ lowerRamificationGroup B G n := by
  rw [lowerRamificationGroup_subtype]
  exact Subgroup.mem_subgroupOf

/-- ★**`i(σ)` は部分群に落としても変わらない**。★`rfl` で閉じる
(`↥H` の作用は `Subgroup.subtype` に沿った制限なので `σ • α` が定義的に同じ)。 -/
theorem ramIndex_subtype {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {H : Subgroup G} (α : B) (τ : ↥H) :
    ramIndex (G := ↥H) α τ = ramIndex (G := G) α (τ : G) := rfl

/-- ★★**段 2 の `h1` を `H = G_1` について供給する** —— `↥G_1` の中で見れば
`(G_1)_1 = ⊤`、すなわち「`↥G_1` は自分自身の第 1 分岐群である」。

`H.subgroupOf H = ⊤`(`Subgroup.subgroupOf_self`)に帰着する。 -/
theorem lowerRamificationGroup_subtype_one_eq_top {B : Type*} [CommRing B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] :
    lowerRamificationGroup B ↥(lowerRamificationGroup B G 1) 1 = ⊤ := by
  rw [lowerRamificationGroup_subtype, Subgroup.subgroupOf_self]

/-- ★★**段 2 の `hne` を `H = G_1` について供給する**(`1 ≤ n` が要る)。

`G_n ≠ G_{n+1}` かつ `1 ≤ n` なら、`↥G_1` の中でも `H_n ≠ H_{n+1}`。

★`1 ≤ n` は落とせない。`n = 0` では `G_0` も `G_1` も `H = G_1` への制限が
`⊤` に潰れるので、左右が等しくなってしまう(ファイル冒頭の注)。 -/
theorem lowerRamificationGroup_subtype_ne {B : Type*} [CommRing B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {n : ℕ} (hn : 1 ≤ n)
    (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) :
    lowerRamificationGroup B ↥(lowerRamificationGroup B G 1) n
      ≠ lowerRamificationGroup B ↥(lowerRamificationGroup B G 1) (n + 1) := by
  rw [lowerRamificationGroup_subtype, lowerRamificationGroup_subtype]
  intro h
  exact hne (subgroupOf_injOn (lowerRamificationGroup_antitone B G hn)
    (lowerRamificationGroup_antitone B G (by omega)) h)

/-! ## §3 ★★★橋の本体 —— `φ_H = φ_{↥H}` -/

/-- ★★★**部分群版 `φ_H` と群版 `φ_{↥H}` は等しい。**

    herbrandPhi α H n = herbrandPhiGroup ↥H α n

★`⊤` の場合(`herbrandPhiGroup_eq_herbrandPhi_top`)は
`Fintype.sum_equiv Subgroup.topEquiv` と `Nat.card_congr` を要したが、
**一般の `H` では `rfl` で閉じる**。どちらも添字型 `↥H` の上の
`phiOf (fun τ => ramIndex α ↑τ)` に展開されるからである(`ramIndex_subtype`)。

★両辺の `[Fintype ↥H]` は同一のインスタンスである(片方は `herbrandPhi` の
`[Fintype ↥H]`、片方は `herbrandPhiGroup` の `[Fintype G]` に `G := ↥H` を
代入したもの)。 -/
theorem herbrandPhi_eq_herbrandPhiGroup_subtype {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    (α : B) (H : Subgroup G) [Fintype ↥H] (n : ℝ) :
    herbrandPhi α H n = herbrandPhiGroup ↥H α n := by
  rw [herbrandPhi_eq_phiOf, herbrandPhiGroup]
  rfl

/-- 存在形の移送。★段 1+2 の結論(`φ_{↥H}(n) ∈ ℤ≥0`)を段 3 が要る形
(`φ_H(n) ∈ ℤ≥0`)に読み替える。 -/
theorem exists_natCast_herbrandPhi_of_herbrandPhiGroup_subtype {B : Type*} [CommRing B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    (α : B) (H : Subgroup G) [Fintype ↥H] (n : ℝ)
    (h : ∃ j : ℕ, herbrandPhiGroup ↥H α n = (j : ℝ)) :
    ∃ k : ℕ, herbrandPhi α H n = (k : ℝ) := by
  obtain ⟨j, hj⟩ := h
  exact ⟨j, by rw [herbrandPhi_eq_herbrandPhiGroup_subtype, hj]⟩

/-! ## §4 ★★★段 3 の `hk` を段 1+2 から供給する -/

/-- ★★★★**段 3 の仮定 `hk` を段 1+2 から供給する**。

原文 (Yoshida08 p.16) が括弧で断っている部分:

> (we know φH(n) ∈Z≥0)

および `n = 0` の場合:

> If n = 0 then φH(0) = 0.

`H = G_1` は `H = H_1` を満たす(`lowerRamificationGroup_subtype_one_eq_top`)ので、
`↥H` に対して**段 2**
(`exists_natCast_herbrandPhiGroup_of_lowerRamificationGroup_one_eq_top`)が
そのまま使える。段 2 の結論は `herbrandPhiGroup ↥H π' n` の形なので、
§3 の橋で `herbrandPhi π' H n` に読み替える。

★仮定は `A`・`B`・`G` についてはすべて段 2 のもの(原典 §6.1 の底の設定)である。
`G_1` について新たに置いた仮定は `[Fintype ↥(lowerRamificationGroup B G 1)]` のみで、
これは段 3 が既に持っているものと同じである。

★`n = 0` は `herbrandPhi_zero` で片付く(段 2 は呼ばない)。 -/
theorem exists_natCast_herbrandPhi_lowerRamificationGroup_one
    {A B : Type*} [CommRing A] [IsDomain A] [IsDiscreteValuationRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] [IsNoetherian A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B]
    [FaithfulSMul G B] [Fintype ↥(lowerRamificationGroup B G 1)]
    (p : ℕ) (hp : p.Prime) [CharP (ResidueField B) p]
    {π' : B} (hπ' : Irreducible π')
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAinj : Function.Injective (algebraMap A B))
    (habel : ∀ x y : G, x * y = y * x)
    {n : ℕ} (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) :
    ∃ k : ℕ, herbrandPhi π' (lowerRamificationGroup B G 1) (n : ℝ) = (k : ℝ) := by
  rcases Nat.eq_zero_or_pos n with rfl | hn
  · exact ⟨0, by
      rw [Nat.cast_zero,
        herbrandPhi_zero ((irreducible_iff_uniformizer π').mp hπ') _]⟩
  · refine exists_natCast_herbrandPhi_of_herbrandPhiGroup_subtype π' _ _ ?_
    exact exists_natCast_herbrandPhiGroup_of_lowerRamificationGroup_one_eq_top
      (A := A) (G := ↥(lowerRamificationGroup B G 1)) p hp hπ' hresA hadj hAinj
      (comm_subtype habel _) lowerRamificationGroup_subtype_one_eq_top
      (lowerRamificationGroup_subtype_ne hn hne)

/-! ## §5 段 3 —— `hk` を仮定から外した形 -/

/-- ★★★**段 3(`G ≠ G_1`)から `hk` を外した形**。

Y15 の `exists_natCast_herbrandPhiGroup_of_tame_quotient` は

    hk : herbrandPhi π (lowerRamificationGroup B G 1) (n : ℝ) = (k : ℝ)

を仮定として持っていた。§4 がそれを供給するので、ここでは仮定から消える。
★**`HasseArf.lean` は 1 行も書き換えていない**(新しい名前で立て直しただけ)。

残る仮定 `hcomp` / `htopC` / `honeC` は順分岐商 `K′^{G_1}/K` の塔に関するもので、
別ノードの担当である。 -/
theorem exists_natCast_herbrandPhiGroup_of_tame_quotient_of_setup
    {A B : Type*} [CommRing A] [IsDomain A] [IsDiscreteValuationRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] [IsNoetherian A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B]
    [FaithfulSMul G B] [Fintype ↥(lowerRamificationGroup B G 1)]
    {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [MulSemiringAction G C]
    (p : ℕ) [hp : Fact p.Prime] [CharP (ResidueField B) p]
    {π' : B} (hπ' : Irreducible π')
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAinj : Function.Injective (algebraMap A B))
    (h0 : lowerRamificationGroup B G 0 = ⊤) (habel : ∀ x y : G, x * y = y * x)
    {n : ℕ} (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) {ϖ : C}
    (hcomp : herbrandPhiGroup G π' (n : ℝ)
      = herbrandPhiGroup G ϖ (herbrandPhi π' (lowerRamificationGroup B G 1) (n : ℝ)))
    (htopC : ∀ σ : G, σ ∈ lowerRamificationGroup B G 1 → ramIndex ϖ σ = ⊤)
    (honeC : ∀ σ : G, σ ∉ lowerRamificationGroup B G 1 → ramIndex ϖ σ = 1) :
    ∃ k' : ℕ, herbrandPhiGroup G π' (n : ℝ) = (k' : ℝ) := by
  obtain ⟨k, hk⟩ := exists_natCast_herbrandPhi_lowerRamificationGroup_one (A := A) p hp.out hπ'
    hresA hadj hAinj habel hne
  exact exists_natCast_herbrandPhiGroup_of_tame_quotient (A := A) (C := C) p hπ' hadj h0 habel
    hne hk hcomp htopC honeC

/-! ## §6 ★★★★★ Hasse-Arf(Theorem 6.11)—— 3 段の合流 -/

/-- ★★★★★**Yoshida 2008 Theorem 6.11 (Hasse-Arf) —— 段 1+2+3 の合流**。

原文 (Yoshida08 p.16):
> Theorem 6.11 (Hasse-Arf). If G is abelian, n ∈ Z[bb]_≥0 and G_n = G_n+1, then φ_G(n) ∈ Z[bb]_≥0.

★逐語の `G_n = G_{n+1}` は `pdftotext` が `≠` の斜線を落とした形である
(原典は `G_n ≠ G_{n+1}`)。

証明は原文どおりの 2 分岐:

* `G = G_1`(`h1 : lowerRamificationGroup B G 1 = ⊤`)—— 段 2
  (`HasseArfStrongInduction.lean`)。★この分岐では `C` も `hC` も使わない。
* `G ≠ G_1` —— 段 3(`HasseArf.lean`)に §4 の `hk` を供給した §5。

★★**`hC` は `G ≠ G_1` のときにだけ要求される**(逸脱の記録 1)。
中身は順分岐商 `K′^{G_1}/K` に対応する DVR `C` と素元 `ϖ` についての 3 条件で、
その構成は別ノードの担当である。★本定理はそこを**仮定として受け取ったまま**であり、
それ以外の仮定はすべて原典 §6.1 の底の設定(`π′` は素元、`𝒪_{K′} = 𝒪_K[π′]`、
`𝒪_K → 𝒪_{K′}` 単射、`K′/K` 完全分岐 `h0`、剰余標数 `p`)と `G` 可換である。 -/
theorem exists_natCast_herbrandPhiGroup_of_abelian
    {A B : Type*} [CommRing A] [IsDomain A] [IsDiscreteValuationRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] [IsNoetherian A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B]
    [FaithfulSMul G B] [Fintype ↥(lowerRamificationGroup B G 1)]
    {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [MulSemiringAction G C]
    (p : ℕ) [hp : Fact p.Prime] [CharP (ResidueField B) p]
    {π' : B} (hπ' : Irreducible π')
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAinj : Function.Injective (algebraMap A B))
    (h0 : lowerRamificationGroup B G 0 = ⊤) (habel : ∀ x y : G, x * y = y * x)
    {n : ℕ} (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) {ϖ : C}
    (hC : lowerRamificationGroup B G 1 ≠ ⊤ →
      herbrandPhiGroup G π' (n : ℝ)
          = herbrandPhiGroup G ϖ (herbrandPhi π' (lowerRamificationGroup B G 1) (n : ℝ))
        ∧ (∀ σ : G, σ ∈ lowerRamificationGroup B G 1 → ramIndex ϖ σ = ⊤)
        ∧ (∀ σ : G, σ ∉ lowerRamificationGroup B G 1 → ramIndex ϖ σ = 1)) :
    ∃ k : ℕ, herbrandPhiGroup G π' (n : ℝ) = (k : ℝ) := by
  by_cases h1 : lowerRamificationGroup B G 1 = ⊤
  · exact exists_natCast_herbrandPhiGroup_of_lowerRamificationGroup_one_eq_top (A := A) p hp.out
      hπ' hresA hadj hAinj habel h1 hne
  · obtain ⟨hcomp, htopC, honeC⟩ := hC h1
    exact exists_natCast_herbrandPhiGroup_of_tame_quotient_of_setup (A := A) (C := C) p hπ'
      hresA hadj hAinj h0 habel hne hcomp htopC honeC

/-! ## §7 余力 —— 段 3 の `htopC` だけは今すぐ供給できる

★段 3 が残している 3 条件のうち `htopC` は**固定環の定義から直ちに出る**。
残る `hcomp`(Lemma 6.10(ii) の 7 つの適合条件)と `honeC`(順分岐商の
第 1 分岐群が自明であること)は別ノードの担当である。 -/

/-- ★**段 3 の `htopC` の供給**(`C = 𝒪_{K′}^H` の場合)。

`H` の元は固定環 `fixedRing B H` の元をすべて固定するので、
その上での `i(σ) = v(σϖ − ϖ)` は `⊤` である。

★★`ramIndex_eq_top_iff` と `mem_fixedRing` を 1 回ずつ使うだけで、
分岐の議論は 1 つも要らない。★`[H.Normal]` は
`MulSemiringAction G ↥(fixedRing B H)`(`fixedRingMulSemiringAction`)を
出すために要る(`lean-idioms.md` #125(iii))。 -/
theorem ramIndex_fixedRing_eq_top_of_mem {B : Type*} [CommRing B] [IsDomain B] {G : Type*}
    [Group G] [MulSemiringAction G B] {H : Subgroup G} [H.Normal]
    [IsDiscreteValuationRing ↥(fixedRing B H)]
    (ϖ : ↥(fixedRing B H)) {σ : G} (hσ : σ ∈ H) :
    ramIndex ϖ σ = ⊤ :=
  (ramIndex_eq_top_iff ϖ σ).2 (Subtype.ext (mem_fixedRing.1 ϖ.2 σ hσ))

/-- ★段 1+2+3 の結論の非負性(原文の `φ_G(n) ∈ Z≥0` の `≥0`)。 -/
theorem herbrandPhiGroup_nonneg_of_abelian
    {A B : Type*} [CommRing A] [IsDomain A] [IsDiscreteValuationRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] [IsNoetherian A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B]
    [FaithfulSMul G B] [Fintype ↥(lowerRamificationGroup B G 1)]
    {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [MulSemiringAction G C]
    (p : ℕ) [hp : Fact p.Prime] [CharP (ResidueField B) p]
    {π' : B} (hπ' : Irreducible π')
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAinj : Function.Injective (algebraMap A B))
    (h0 : lowerRamificationGroup B G 0 = ⊤) (habel : ∀ x y : G, x * y = y * x)
    {n : ℕ} (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) {ϖ : C}
    (hC : lowerRamificationGroup B G 1 ≠ ⊤ →
      herbrandPhiGroup G π' (n : ℝ)
          = herbrandPhiGroup G ϖ (herbrandPhi π' (lowerRamificationGroup B G 1) (n : ℝ))
        ∧ (∀ σ : G, σ ∈ lowerRamificationGroup B G 1 → ramIndex ϖ σ = ⊤)
        ∧ (∀ σ : G, σ ∉ lowerRamificationGroup B G 1 → ramIndex ϖ σ = 1)) :
    0 ≤ herbrandPhiGroup G π' (n : ℝ) :=
  nonneg_of_exists_natCast (exists_natCast_herbrandPhiGroup_of_abelian (A := A) (C := C) p hπ'
    hresA hadj hAinj h0 habel hne hC)

end ABC3.Found.PGC
