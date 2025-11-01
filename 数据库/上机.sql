-------------9-------------
SELECT Name, Salary from Employees order by Salary DESC;


------------10--------------
SELECT DISTINCT
    author_id as id
from Views
where
    author_id = viewer_id
order by id;

------------11--------------
SELECT name
from Employee
where
    id in (
        SELECT managerId
        from Employee
        group by
            managerId
        having
            count(id) >= 5
    );

------------12--------------
SELECT s.user_id, ROUND(
        count(
            CASE
                WHEN c.action = 'confirmed' THEN 1
            END
        ) * 1.0 / count(*), 2
    ) AS confirmation_rate
from Signups s
    LEFT JOIN Confirmations c ON s.user_id = c.user_id
group by
    s.user_id;

------------13--------------
SELECT
    DATE_FORMAT(trans_date, '%Y-%m') as month,
    country,
    count(*) as trans_count,
    count(
        CASE
            WHEN state = 'approved' THEN 1
        END
    ) as approved_count,
    sum(amount) as trans_total_amount,
    sum(
        CASE
            WHEN state = 'approved' THEN amount
            else 0
        END
    ) as approved_total_amount
from Transactions
group by
    month,
    country;

------------14--------------
SELECT ROUND(
        SUM(
            CASE
                WHEN d.order_date = d.customer_pref_delivery_date THEN 1
                ELSE 0
            END
        ) * 100.0 / COUNT(DISTINCT d.customer_id), 2
    ) AS immediate_percentage
FROM Delivery d
WHERE (d.customer_id, d.order_date) IN (
        SELECT customer_id, MIN(order_date)
        FROM Delivery
        GROUP BY
            customer_id
    );

------------15--------------
SELECT activity_date as date, count(DISTINCT user_id) as active_users
from Activity
where
    activity_date <= '2024-07-27'
    and DATEDIFF(activity_date, '2024-07-27') > -30
group by
    date;

------------16--------------
(
    select num as max_single_number
    from MyNumbers
    group by
        num
    having
        count(*) = 1
    order by num DESC
    limit 1
)
union
(
    SELECT null
)
limit 1;

------------17--------------
SELECT user_id, count(*) as friend_count
from (
        SELECT
            requester_id as user_id, accepter_id as f_id
        from RequestAccepted
        union ALL
        SELECT
            accepter_id as user_id, requester_id as f_id
        from RequestAccepted
    ) AS temp
group by
    user_id
order by friend_count DESC
limit 1;

------------18--------------
SELECT person_name
from (
        SELECT person_name, sum(weight) over (
                order by turn
            ) as s
        from Queue
    ) as temp
where
    s <= 1000
order by s DESC
limit 1;

------------19--------------
SELECT u.user_id as buyer_id, u.join_date, count(o.order_id) as orders_in_2019
from (
        SELECT user_id, join_date
        from Users
    ) u
    LEFT join (
        SELECT buyer_id, order_id
        from Orders
        where
            year(order_date) = 2019
    ) o on u.user_id = o.buyer_id
group by
    u.user_id;